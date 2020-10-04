import numpy as np 
import unittest
from geograd import geograd_parallel as g
from geograd_complex import geograd_parallel as gcs
from stl import mesh
import os
h = 1e-15
from openmdao.utils.assert_utils import assert_rel_error
from numpy.testing import assert_almost_equal
import warnings 
from mpi4py import MPI

def custom_assert(self, truth, approx, base_tol=1e-7):
    if np.abs(truth) > 0.1:
        assert_rel_error(self, truth, approx, tolerance=base_tol)
    elif np.abs(truth) > 1e-4:
        assert_rel_error(self, truth, approx, tolerance=base_tol*10)
    elif np.abs(truth) > 5e-9:
        assert_rel_error(self, truth, approx, tolerance=base_tol*100)
    else:
        assert_almost_equal(truth, approx, decimal=7)

def test_derivatives_CS(A1, B1, C1, A2, B2, C2, rho, testcase, indices_1=None, indices_2=None, method='cs'):
        maxdim = np.max(np.maximum(np.maximum(A2.max(axis=1), B2.max(axis=1)), C2.max(axis=1)) - np.minimum(np.minimum(A2.min(axis=1), B2.min(axis=1)), C2.min(axis=1)))
        result = g.compute(A1, B1, C1, A2, B2, C2, 1.0, rho, maxdim)
        result2 = g.compute_derivs(A1, B1, C1, A2, B2, C2, result[2], rho, maxdim)
        ks_base = result2[0]
        A1_grad = result2[5]
        B1_grad = result2[6]
        C1_grad = result2[7]
        A2_grad = result2[8]
        B2_grad = result2[9]
        C2_grad = result2[10]

        # because of the way STL files are stored, mesh vertices appear multiple times in the ABC matrices.
        # they need to be manipulated at the same time or the minimum distance gradients won't be correct

        cseps = 1.0e-10
        fdeps = 5.0e-3
        
        if indices_1 is None:
            indices_1 = [0]
        if indices_2 is None:
            indices_2 = [0]
        mesh1_contrib_indices = np.array(indices_1,dtype=np.int_)
        mesh2_contrib_indices = np.array(indices_2,dtype=np.int_)

        A1_pts = np.take(A1, mesh1_contrib_indices, axis=1)
        B1_pts = np.take(B1, mesh1_contrib_indices, axis=1)
        C1_pts = np.take(C1, mesh1_contrib_indices, axis=1)
        p1_to_perturb = np.unique(np.hstack([A1_pts, B1_pts, C1_pts]), axis=1)

        tot_counter = 0
        nonzero_grad_ks = False

        #first perturb A1 B1 C1 mesh points
        for i in range(p1_to_perturb.shape[1]):
            # search each input for a particular point and get the indices where that point appears
            vec_to_search_for = p1_to_perturb[:,i].reshape((3,1))
            # print('Search vec: '+str(vec_to_search_for))
            A1_indices = np.argwhere(np.all((A1-vec_to_search_for)==0, axis=0))
            B1_indices = np.argwhere(np.all((B1-vec_to_search_for)==0, axis=0))
            C1_indices = np.argwhere(np.all((C1-vec_to_search_for)==0, axis=0))

        #     # 4) do a complex step increment to each x, y, and z component for each point at all locations in p and q where it appears
        #     # this alters the mesh in a self-consistent fashion same as pyGeo would
            if method=='cs':
                for k in range(3):
                    exact_grad_ks = 0
                    A1_alt = A1.copy().astype(np.complex64)
                    for j in A1_indices:
                        jind = j[0]
                        A1_alt[k,jind] = A1_alt[k,jind] + cseps*1.0j
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += A1_grad[k,jind]
                    B1_alt = B1.copy().astype(np.complex64)
                    for j in B1_indices:
                        jind = j[0]
                        B1_alt[k,jind] = B1_alt[k,jind] + cseps*1.0j
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += B1_grad[k,jind]
                    C1_alt = C1.copy().astype(np.complex64)
                    for j in C1_indices:
                        jind = j[0]
                        C1_alt[k,jind] = C1_alt[k,jind] + cseps*1.0j
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += C1_grad[k,jind]
                    resultcs = gcs.compute(A1_alt, B1_alt, C1_alt, A2, B2, C2, 1.0, rho, maxdim)
                    resultcs2 = gcs.compute(A1_alt, B1_alt, C1_alt, A2, B2, C2, resultcs[2], rho, maxdim)
                    ks_cs = resultcs2[0]
                    gradcs_ks = np.imag((ks_cs - ks_base)) / cseps

                    custom_assert(testcase, gradcs_ks, exact_grad_ks)
                    # warnings.warn('Surf mesh CS: '+str(gradcs_ks)+' Exact: '+str(exact_grad_ks))
                    tot_counter += 1
                    if np.abs(gradcs_ks) > 1e-3:
                        nonzero_grad_ks = True

            else:
                for k in range(3):
                    exact_grad_ks = 0
                    A1_alt = A1.copy()
                    for j in A1_indices:
                        jind = j[0]
                        A1_alt[k,jind] = A1_alt[k,jind] + fdeps
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += A1_grad[k,jind]
                    B1_alt = B1.copy()
                    for j in B1_indices:
                        jind = j[0]
                        B1_alt[k,jind] = B1_alt[k,jind] + fdeps
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += B1_grad[k,jind]
                    C1_alt = C1.copy()
                    for j in C1_indices:
                        jind = j[0]
                        C1_alt[k,jind] = C1_alt[k,jind] + fdeps
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += C1_grad[k,jind]

                    resultfd = g.compute(A1_alt, B1_alt, C1_alt, A2, B2, C2, 1.0, rho, maxdim)
                    resultfd2 = g.compute(A1_alt, B1_alt, C1_alt, A2, B2, C2, resultfd[2], rho, maxdim)
                    ks_fd = resultfd2[0]
                    gradfd_ks = (ks_fd - ks_base) / fdeps

                    custom_assert(testcase, gradfd_ks, exact_grad_ks)


        testcase.assertTrue(tot_counter > 0, msg='If tot_counter remains zero there is some issue finding close points using the numpy slicing and searching')
        testcase.assertTrue(nonzero_grad_ks, msg='Check to make sure at least one gradient checked is actually nonzero')

        A2_pts = np.take(A2, mesh2_contrib_indices, axis=1)
        B2_pts = np.take(B2, mesh2_contrib_indices, axis=1)
        C2_pts = np.take(C2, mesh2_contrib_indices, axis=1)
        p2_to_perturb = np.unique(np.hstack([A2_pts, B2_pts, C2_pts]), axis=1)

        tot_counter = 0
        nonzero_grad_ks = False

        #next perturb surface mesh points
        for i in range(p2_to_perturb.shape[1]):
            # search each input for a particular point and get the indices where that point appears
            vec_to_search_for = p2_to_perturb[:,i].reshape((3,1))
            # print('Search vec: '+str(vec_to_search_for))
            A2_indices = np.argwhere(np.all((A2-vec_to_search_for)==0, axis=0))
            B2_indices = np.argwhere(np.all((B2-vec_to_search_for)==0, axis=0))
            C2_indices = np.argwhere(np.all((C2-vec_to_search_for)==0, axis=0))

        #     # 4) do a complex step increment to each x, y, and z component for each point at all locations in p and q where it appears
        #     # this alters the mesh in a self-consistent fashion same as pyGeo would
            if method=='cs':
                for k in range(3):
                    exact_grad_ks = 0
                    A2_alt = A2.copy().astype(np.complex64)
                    for j in A2_indices:
                        jind = j[0]
                        A2_alt[k,jind] = A2_alt[k,jind] + cseps*1.0j
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += A2_grad[k,jind]
                    B2_alt = B2.copy().astype(np.complex64)
                    for j in B2_indices:
                        jind = j[0]
                        B2_alt[k,jind] = B2_alt[k,jind] + cseps*1.0j
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += B2_grad[k,jind]
                    C2_alt = C2.copy().astype(np.complex64)
                    for j in C2_indices:
                        jind = j[0]
                        C2_alt[k,jind] = C2_alt[k,jind] + cseps*1.0j
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += C2_grad[k,jind]

                    resultcs = gcs.compute(A1, B1, C1, A2_alt, B2_alt, C2_alt, 1.0, rho, maxdim)
                    resultcs2 = gcs.compute(A1, B1, C1, A2_alt, B2_alt, C2_alt, resultcs[2], rho, maxdim)
                    ks_cs = resultcs2[0]
                    gradcs_ks = np.imag((ks_cs - ks_base)) / cseps

                    custom_assert(testcase, gradcs_ks, exact_grad_ks)
                    # warnings.warn('Surf mesh CS: '+str(gradcs_ks)+' Exact: '+str(exact_grad_ks))
                    tot_counter += 1
                    if np.abs(gradcs_ks) > 1e-3:
                        nonzero_grad_ks = True

            else:
                for k in range(3):
                    exact_grad_ks = 0
                    A2_alt = A2.copy()
                    for j in A2_indices:
                        jind = j[0]
                        A2_alt[k,jind] = A2_alt[k,jind] + fdeps
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += A2_grad[k,jind]
                    B2_alt = B2.copy()
                    for j in B2_indices:
                        jind = j[0]
                        B2_alt[k,jind] = B2_alt[k,jind] + fdeps
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += B2_grad[k,jind]
                    C2_alt = C2.copy()
                    for j in C2_indices:
                        jind = j[0]
                        C2_alt[k,jind] = C2_alt[k,jind] + fdeps
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += C2_grad[k,jind]

                    resultfd = g.compute(A1, B1, C1, A2_alt, B2_alt, C2_alt, 1.0, rho, maxdim)
                    resultfd2 = g.compute(A1, B1, C1, A2_alt, B2_alt, C2_alt, resultfd[2], rho, maxdim)
                    ks_fd = resultfd2[0]
                    gradfd_ks = (ks_fd - ks_base) / fdeps

                    custom_assert(testcase, gradfd_ks, exact_grad_ks)

        testcase.assertTrue(tot_counter > 0, msg='If tot_counter remains zero there is some issue finding close points using the numpy slicing and searching')
        testcase.assertTrue(nonzero_grad_ks, msg='Check to make sure at least one gradient checked is actually nonzero')

class TestIntersectTrivial(unittest.TestCase):
    def setUp(self):
        self.a1 = np.repeat(np.array([0.0, 0.0, 0.0]),2).reshape(3,2)
        self.b1 = np.repeat(np.array([1.0, 0.0, 0.0]),2).reshape(3,2)
        self.c1 = np.repeat(np.array([0.0, 1.0, 0.0]),2).reshape(3,2)

    def test_intersect_1_in_2(self):
        a2 = np.repeat(np.array([-1.0, 0.5, -0.1]),2).reshape(3,2)
        b2 = np.repeat(np.array([1.0, 0.5, -0.1]),2).reshape(3,2)
        c2 = np.repeat(np.array([0.5, 0.5, 5.0]),2).reshape(3,2)
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0, 10, 1000.0)
        self.assertAlmostEqual(result[1], 0.5*4, 10)

    def test_intersect_2_in_1(self):
        a2 = np.repeat(np.array([0.75, 0.5, -2.0]),2).reshape(3,2)
        b2 = np.repeat(np.array([-0.25, 0.5, -2.0]),2).reshape(3,2)
        c2 = np.repeat(np.array([0.25, 0.5, 1.0]),2).reshape(3,2)
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0, 10, 1000.0)
        self.assertAlmostEqual(result[1], 1/3*4, 10)

    def test_intersect_right_edge(self):
        a2 = np.repeat(np.array([1.0, 0.5, -2.0]),2).reshape(3,2)
        b2 = np.repeat(np.array([0.0, 0.5, -2.0]),2).reshape(3,2)
        c2 = np.repeat(np.array([0.5, 0.5, 1.0]),2).reshape(3,2)
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0, 10, 1000.0)
        self.assertAlmostEqual(result[1], 1/6*4, 10)

    def test_intersect_left_edge(self):
        a2 = np.repeat(np.array([0.5, 0.5, -2.0]),2).reshape(3,2)
        b2 = np.repeat(np.array([-0.5, 0.5, -2.0]),2).reshape(3,2)
        c2 = np.repeat(np.array([0.0, 0.5, 1.0]),2).reshape(3,2)
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0, 10, 1000.0)
        self.assertAlmostEqual(result[1], 1/6*4, 10)

class MinDistSTLTestCase1(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.getcwd()+'/tests'):
            test_data_path = os.getcwd()+'/tests'
        else:
            raise IOError

        # the 'surface'  mesh is a the same blob file offset by 3 units in the y direction
        surface_mesh = mesh.Mesh.from_file(test_data_path+'/blob_offset_y3.stl').vectors
        object_mesh = mesh.Mesh.from_file(test_data_path+'/blob.stl').vectors
        self.smSize = surface_mesh[:,0,:].shape[0]
        self.smp0 = surface_mesh[:,0,:].transpose()
        self.smp1 = surface_mesh[:,1,:].transpose()
        self.smp2 = surface_mesh[:,2,:].transpose()

        self.omSize = object_mesh[:,0,:].shape[0]
        self.objp0 = object_mesh[:,0,:].transpose()
        self.objp1 = object_mesh[:,1,:].transpose()
        self.objp2 = object_mesh[:,2,:].transpose()
        self.maxdim = np.max(np.maximum(np.maximum(self.objp0.max(axis=1), self.objp1.max(axis=1)), self.objp2.max(axis=1)) - np.minimum(np.minimum(self.objp0.min(axis=1), self.objp1.min(axis=1)), self.objp2.min(axis=1)))


    def test_values(self):
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 10, self.maxdim)
        self.assertAlmostEqual(2.33166, result[2], 4)
        self.assertAlmostEqual(0.0, result[1], 10)
        result2 = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300, self.maxdim)
        self.assertAlmostEqual(-2.3197215104930256, result2[0], 8)
    
    def test_derivs(self):
        # closepoint is at 40, 47
        test_derivatives_CS(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 300, self, [40, 10, 91], [47, 5, 90])



class MinDistSTLTestCase2(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.getcwd()+'/tests'):
            test_data_path = os.getcwd()+'/tests'
        else:
            raise IOError

        # the 'surface'  mesh is a the same blob file offset by 3 units in the y direction
        surface_mesh = mesh.Mesh.from_file(test_data_path+'/blob_offset_y3_line_line.stl').vectors
        object_mesh = mesh.Mesh.from_file(test_data_path+'/blob_line_line.stl').vectors
        self.smSize = surface_mesh[:,0,:].shape[0]
        self.smp0 = surface_mesh[:,0,:].transpose()
        self.smp1 = surface_mesh[:,1,:].transpose()
        self.smp2 = surface_mesh[:,2,:].transpose()

        self.omSize = object_mesh[:,0,:].shape[0]
        self.objp0 = object_mesh[:,0,:].transpose()
        self.objp1 = object_mesh[:,1,:].transpose()
        self.objp2 = object_mesh[:,2,:].transpose()
        self.maxdim = np.max(np.maximum(np.maximum(self.objp0.max(axis=1), self.objp1.max(axis=1)), self.objp2.max(axis=1)) - np.minimum(np.minimum(self.objp0.min(axis=1), self.objp1.min(axis=1)), self.objp2.min(axis=1)))

    def test_values(self):
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 10, self.maxdim)
        self.assertAlmostEqual(0.0178387, result[2], 6)
        self.assertAlmostEqual(0.0, result[1], 10)
        result2 = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300, self.maxdim)
        self.assertAlmostEqual(-0.013207541564420895, result2[0], 8)

    def test_derivs(self):
        test_derivatives_CS(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 300, self, [61, 10, 87], [51, 10, 12])

class MinDistSTLTestCase3(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.getcwd()+'/tests'):
            test_data_path = os.getcwd()+'/tests'
        else:
            raise IOError

        # the 'surface'  mesh is a the same blob file offset by 3 units in the y direction
        surface_mesh = mesh.Mesh.from_file(test_data_path+'/blob_offset_y3_point_tri.stl').vectors
        object_mesh = mesh.Mesh.from_file(test_data_path+'/blob_point_tri.stl').vectors
        self.smSize = surface_mesh[:,0,:].shape[0]
        self.smp0 = surface_mesh[:,0,:].transpose()
        self.smp1 = surface_mesh[:,1,:].transpose()
        self.smp2 = surface_mesh[:,2,:].transpose()

        self.omSize = object_mesh[:,0,:].shape[0]
        self.objp0 = object_mesh[:,0,:].transpose()
        self.objp1 = object_mesh[:,1,:].transpose()
        self.objp2 = object_mesh[:,2,:].transpose()
        self.maxdim = np.max(np.maximum(np.maximum(self.objp0.max(axis=1), self.objp1.max(axis=1)), self.objp2.max(axis=1)) - np.minimum(np.minimum(self.objp0.min(axis=1), self.objp1.min(axis=1)), self.objp2.min(axis=1)))

    def test_values(self):
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 10, self.maxdim)
        self.assertAlmostEqual(0.02212236, result[2], 6)
        self.assertAlmostEqual(0.0, result[1], 10)
        result2 = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300, self.maxdim)
        self.assertAlmostEqual(-0.015189136030326717, result2[0], 8)
    
    def test_derivs(self):
        test_derivatives_CS(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 300, self, [104, 87], [52, 2])

if __name__ == '__main__':
    unittest.main()
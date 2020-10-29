import numpy as np 
import unittest
from geograd import geograd_serial as g
from geograd import geograd_serial_complex as gcs
from stl import mesh
import os
h = 1e-15
from openmdao.utils.assert_utils import assert_near_equal
from numpy.testing import assert_almost_equal
import warnings 

def custom_assert(self, truth, approx, base_tol=1e-7):
    if np.abs(truth) > 0.1:
        assert_near_equal(truth, approx, tolerance=base_tol)
    elif np.abs(truth) > 1e-4:
        assert_near_equal(truth, approx, tolerance=base_tol*10)
    elif np.abs(truth) > 5e-9:
        assert_near_equal(truth, approx, tolerance=base_tol*100)
    else:
        assert_almost_equal(truth, approx, decimal=7)

def helper_test_derivatives_translate_objects_random(testcase, objp0, objp1, objp2, smp0, smp1, smp2, n):
        result = g.compute(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        perim = result[1]

        # test derivs
        result2 = g.compute_derivs(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        A1_grad = result2[9]
        B1_grad = result2[10]
        C1_grad = result2[11]
        partial_sum_1  = np.sum(A1_grad, axis=1) + np.sum(B1_grad, axis=1) + np.sum(C1_grad, axis=1)

        A2_grad = result2[12]
        B2_grad = result2[13]
        C2_grad = result2[14]
        partial_sum_2 = np.sum(A2_grad, axis=1) + np.sum(B2_grad, axis=1) + np.sum(C2_grad, axis=1)
        cseps = 1e-15
        max_abs_der = 0.0
        for i in range(n):
            # test translating derivatives for cube1
            offsetdir1 = np.random.uniform(-1, 1, size=(3,1))
            offsetdir1 = offsetdir1 / np.linalg.norm(offsetdir1)
            objp0cs = objp0 + offsetdir1*cseps*1j
            objp1cs = objp1 + offsetdir1*cseps*1j
            objp2cs = objp2 + offsetdir1*cseps*1j
            rescs = gcs.compute(objp0cs, objp1cs, objp2cs, smp0, smp1, smp2, 0.001, 10)
            gradcs = np.imag((rescs[1] - perim)) / cseps
            gradexact = np.sum(offsetdir1*partial_sum_1.reshape(3,1))

            assert_almost_equal(gradexact, gradcs, decimal=10)
            if np.abs(gradcs) > max_abs_der:
                max_abs_der = np.abs(gradcs)

            # test translating derivs for cube2
            offsetdir2 = np.random.uniform(-1, 1, size=(3,1))
            smp0cs = smp0 + offsetdir2*cseps*1j
            smp1cs = smp1 + offsetdir2*cseps*1j
            smp2cs = smp2 + offsetdir2*cseps*1j
            rescs = gcs.compute(objp0, objp1, objp2, smp0cs, smp1cs, smp2cs, 0.001, 10)
            gradcs = np.imag((rescs[1] - perim)) / cseps
            gradexact = np.sum(offsetdir2*partial_sum_2.reshape(3,1))

            assert_almost_equal(gradexact, gradcs, decimal=10)
            if np.abs(gradcs) > max_abs_der:
                max_abs_der = np.abs(gradcs)
        testcase.assertGreater(max_abs_der, 1e-8)

def helper_test_derivatives_translate_object_fd(testcase, objp0, objp1, objp2, smp0, smp1, smp2, direction, value, test_first=True):
        result = g.compute(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        perim = result[1]

        # test derivs
        result2 = g.compute_derivs(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        A1_grad = result2[9]
        B1_grad = result2[10]
        C1_grad = result2[11]
        partial_sum_1  = np.sum(A1_grad, axis=1) + np.sum(B1_grad, axis=1) + np.sum(C1_grad, axis=1)

        A2_grad = result2[12]
        B2_grad = result2[13]
        C2_grad = result2[14]
        partial_sum_2 = np.sum(A2_grad, axis=1) + np.sum(B2_grad, axis=1) + np.sum(C2_grad, axis=1)
        fdeps = 1e-5

        # test translating derivatives for object 1
        if test_first:
            offsetdir1 = direction.copy().reshape((3,1))
            objp0fd = objp0 + offsetdir1*fdeps
            objp1fd = objp1 + offsetdir1*fdeps
            objp2fd = objp2 + offsetdir1*fdeps
            resfd = g.compute(objp0fd, objp1fd, objp2fd, smp0, smp1, smp2, 0.001, 10)
            gradfd = (resfd[1] - perim) / fdeps
            gradexact = np.sum(offsetdir1*partial_sum_1.reshape(3,1))

            assert_almost_equal(gradexact, gradfd, decimal=10)
            assert_almost_equal(gradexact, value)

        # test translating derivs for object 2
        offsetdir2 = - direction.copy().reshape((3,1))
        smp0fd = smp0 + offsetdir2*fdeps
        smp1fd = smp1 + offsetdir2*fdeps
        smp2fd = smp2 + offsetdir2*fdeps
        resfd = g.compute(objp0, objp1, objp2, smp0fd, smp1fd, smp2fd, 0.001, 10)
        gradfd = (resfd[1] - perim) / fdeps
        gradexact = np.sum(offsetdir2*partial_sum_2.reshape(3,1))

        assert_almost_equal(gradexact, gradfd, decimal=10)
        assert_almost_equal(gradexact, value)


def helper_test_derivatives_translate_object_given(testcase, objp0, objp1, objp2, smp0, smp1, smp2, direction, value, test_first=True, tol=1e-7):
        result = g.compute(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        perim = result[1]

        # test derivs
        result2 = g.compute_derivs(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        A1_grad = result2[9]
        B1_grad = result2[10]
        C1_grad = result2[11]
        partial_sum_1  = np.sum(A1_grad, axis=1) + np.sum(B1_grad, axis=1) + np.sum(C1_grad, axis=1)

        A2_grad = result2[12]
        B2_grad = result2[13]
        C2_grad = result2[14]
        partial_sum_2 = np.sum(A2_grad, axis=1) + np.sum(B2_grad, axis=1) + np.sum(C2_grad, axis=1)
        cseps = 1e-15

        # test translating derivatives for object 1
        if test_first:
            offsetdir1 = direction.copy().reshape((3,1))
            objp0cs = objp0 + offsetdir1*cseps*1j
            objp1cs = objp1 + offsetdir1*cseps*1j
            objp2cs = objp2 + offsetdir1*cseps*1j
            rescs = gcs.compute(objp0cs, objp1cs, objp2cs, smp0, smp1, smp2, 0.001, 10)
            gradcs = np.imag((rescs[1] - perim)) / cseps
            gradexact = np.sum(offsetdir1*partial_sum_1.reshape(3,1))

            assert_almost_equal(gradexact, gradcs, decimal=10)
            custom_assert(testcase, gradexact, value, tol)

        # test translating derivs for object 2
        offsetdir2 = - direction.copy().reshape((3,1))
        smp0cs = smp0 + offsetdir2*cseps*1j
        smp1cs = smp1 + offsetdir2*cseps*1j
        smp2cs = smp2 + offsetdir2*cseps*1j
        rescs = gcs.compute(objp0, objp1, objp2, smp0cs, smp1cs, smp2cs, 0.001, 10)
        gradcs = np.imag((rescs[1] - perim)) / cseps
        gradexact = np.sum(offsetdir2*partial_sum_2.reshape(3,1))

        assert_almost_equal(gradexact, gradcs, decimal=10)
        custom_assert(testcase, gradexact, value, tol)

def helper_test_derivatives_CS(A1, B1, C1, A2, B2, C2, rho, testcase, indices_1=None, indices_2=None, method='cs'):
        result = g.compute(A1, B1, C1, A2, B2, C2, 1.0, rho)
        result2 = g.compute_derivs(A1, B1, C1, A2, B2, C2, result[2], rho)
        ks_base = result2[0]
        A1_grad = result2[3]
        B1_grad = result2[4]
        C1_grad = result2[5]
        A2_grad = result2[6]
        B2_grad = result2[7]
        C2_grad = result2[8]

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

                    resultcs = gcs.compute(A1_alt, B1_alt, C1_alt, A2, B2, C2, 1.0, rho)
                    resultcs2 = gcs.compute(A1_alt, B1_alt, C1_alt, A2, B2, C2, resultcs[2], rho)
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

                    resultfd = g.compute(A1_alt, B1_alt, C1_alt, A2, B2, C2, 1.0, rho)
                    resultfd2 = g.compute(A1_alt, B1_alt, C1_alt, A2, B2, C2, resultfd[2], rho)
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

                    resultcs = gcs.compute(A1, B1, C1, A2_alt, B2_alt, C2_alt, 1.0, rho)
                    resultcs2 = gcs.compute(A1, B1, C1, A2_alt, B2_alt, C2_alt, resultcs[2], rho)
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

                    resultfd = g.compute(A1, B1, C1, A2_alt, B2_alt, C2_alt, 1.0, rho)
                    resultfd2 = g.compute(A1, B1, C1, A2_alt, B2_alt, C2_alt, resultfd[2], rho)
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
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0, 10)
        self.assertAlmostEqual(result[1], 0.5*4, 10)

    def test_intersect_2_in_1(self):
        a2 = np.repeat(np.array([0.75, 0.5, -2.0]),2).reshape(3,2)
        b2 = np.repeat(np.array([-0.25, 0.5, -2.0]),2).reshape(3,2)
        c2 = np.repeat(np.array([0.25, 0.5, 1.0]),2).reshape(3,2)
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0, 10)
        self.assertAlmostEqual(result[1], 1/3*4, 10)

    def test_intersect_right_edge(self):
        a2 = np.repeat(np.array([1.0, 0.5, -2.0]),2).reshape(3,2)
        b2 = np.repeat(np.array([0.0, 0.5, -2.0]),2).reshape(3,2)
        c2 = np.repeat(np.array([0.5, 0.5, 1.0]),2).reshape(3,2)
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0, 10)
        self.assertAlmostEqual(result[1], 1/6*4, 10)

    def test_intersect_left_edge(self):
        a2 = np.repeat(np.array([0.5, 0.5, -2.0]),2).reshape(3,2)
        b2 = np.repeat(np.array([-0.5, 0.5, -2.0]),2).reshape(3,2)
        c2 = np.repeat(np.array([0.0, 0.5, 1.0]),2).reshape(3,2)
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0, 10)
        self.assertAlmostEqual(result[1], 1/6*4, 10)

class MinDistSTLTestCase1(unittest.TestCase):
    def setUp(self):
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        test_data_path = self.base_path + r'/inputFiles'

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

    def test_values(self):
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 10)
        self.assertAlmostEqual(2.33166, result[2], 4)
        self.assertAlmostEqual(0.0, result[1], 10)
        result2 = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300)
        self.assertAlmostEqual(-2.3197215104930256, result2[0], 8)
    
    def test_derivs(self):
        # closepoint is at 40, 47
        helper_test_derivatives_CS(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 300, self, [40, 10, 91], [47, 5, 90])



class MinDistSTLTestCase2(unittest.TestCase):
    def setUp(self):
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        test_data_path = self.base_path + r'/inputFiles'

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

    def test_values(self):
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 10)
        self.assertAlmostEqual(0.0178387, result[2], 6)
        self.assertAlmostEqual(0.0, result[1], 10)
        result2 = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300)
        self.assertAlmostEqual(-0.013207541564420895, result2[0], 8)

    def test_derivs(self):
        helper_test_derivatives_CS(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 300, self, [61, 10, 87], [51, 10, 12])

class MinDistSTLTestCase3(unittest.TestCase):
    def setUp(self):
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        test_data_path = self.base_path + r'/inputFiles'

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

    def test_values(self):
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 10)
        self.assertAlmostEqual(0.02212236, result[2], 6)
        self.assertAlmostEqual(0.0, result[1], 10)
        result2 = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300)
        self.assertAlmostEqual(-0.015189136030326717, result2[0], 8)
    
    def test_derivs(self):
        helper_test_derivatives_CS(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 300, self, [104, 87], [52, 2])

def generate_plane(vector1, vector2, start_point):
    vertex_1 = start_point
    vertex_2 = start_point + vector1
    vertex_3 = start_point + vector2
    vertex_4 = start_point + vector1 + vector2
    p0 = np.vstack([vertex_1, vertex_1]).transpose()
    p1 = np.vstack([vertex_2, vertex_3]).transpose()
    p2 = np.vstack([vertex_4, vertex_4]).transpose()
    return p0, p1, p2

class BisectSphereTestCase(unittest.TestCase):
    def test_bisect_sphere(self):
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        test_data_path = self.base_path + r'/inputFiles'
        surface_mesh_base = mesh.Mesh.from_file(test_data_path+'/25mmsphere.stl').vectors
        smSize = surface_mesh_base[:,0,:].shape[0]

        surface_mesh = surface_mesh_base.copy()
        smp0 = surface_mesh[:,0,:].transpose()
        smp1 = surface_mesh[:,1,:].transpose()
        smp2 = surface_mesh[:,2,:].transpose()

        objp0, objp1, objp2 = generate_plane(np.array([80., 0., 0.]), np.array([0., 80., 0.]), np.array([-40., -40., 0.00001]))
        result = g.compute(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)

        custom_assert(self, result[1], np.pi*25., base_tol=5e-4)
        result2 = g.compute_derivs(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        custom_assert(self, result2[1], np.pi*25., base_tol=5e-4)

        # pick a new offset
        radius = 25/2
        offset = radius * 0.5
        
        objp0, objp1, objp2 = generate_plane(np.array([80., 0., 0.]), np.array([0., 80., 0.]), np.array([-40., -40., 1e-4+offset]))
        result = g.compute(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        custom_assert(self, result[1], 2*np.pi*radius*np.sqrt(1-(offset/radius)**2), base_tol=5e-4)
        deriv = 2*np.pi*radius*(1-(offset/radius)**2)**(-1/2)*(-offset)/radius**2
        # the error is discretization error in the sphere, not necessarily error in the computed derivatives of the STL
        helper_test_derivatives_translate_object_given(self, objp0, objp1, objp2, smp0, smp1, smp2, np.array([[0., 0., 1.]]), deriv, test_first=True, tol=1e-2)

class BisectPlaneTestCase(unittest.TestCase):
    def test_bisect_planes(self):
        
        objp0, objp1, objp2 = generate_plane(np.array([0., 0., 80.]), np.array([0., 80., 0.]), np.array([0.0, -40., -40]))
        smp0, smp1, smp2 = generate_plane(np.array([0., 1.5, 0.0]), np.array([20., 0., 1.]), np.array([-10., 0.0, 0.0]))
        result = g.compute(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        custom_assert(self, result[1], 1.5, base_tol=1e-3)

        smp0, smp1, smp2 = generate_plane(np.array([0., 100, 0.0]), np.array([20., 0., 1.]), np.array([-10., -50.0, 0.0]))
        result = g.compute(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        custom_assert(self, result[1], 80., base_tol=1e-3)

class BisectCubeTestCase(unittest.TestCase):
    def test_plane_cube(self):
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        test_data_path = self.base_path + r'/inputFiles'
        surface_mesh_base = mesh.Mesh.from_file(test_data_path+'/bigcube.stl').vectors
        smSize = surface_mesh_base[:,0,:].shape[0]

        surface_mesh = surface_mesh_base.copy()
        smp0 = surface_mesh[:,0,:].transpose()
        smp1 = surface_mesh[:,1,:].transpose()
        smp2 = surface_mesh[:,2,:].transpose()

        objp0, objp1, objp2 = generate_plane(np.array([0., 0., 60.]), np.array([0., 80., 0.]), np.array([0.01, -40., -40]))
        result = g.compute(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        result2 = g.compute_derivs(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        custom_assert(self, result[1], 16., base_tol=1e-7)
        helper_test_derivatives_translate_object_given(self, objp0, objp1, objp2, smp0, smp1, smp2, np.array([[1., 0., 0.]]), 0.0)

class OffsetCubesTestCase(unittest.TestCase):
    def test_two_cubes(self):
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        test_data_path = self.base_path + r'/inputFiles'
        surface_mesh_base = mesh.Mesh.from_file(test_data_path+'/bigcube.stl').vectors
        object_mesh_base = mesh.Mesh.from_file(test_data_path+'/littlecube.stl').vectors
        smSize = surface_mesh_base[:,0,:].shape[0]
        omSize = object_mesh_base[:,0,:].shape[0]

        surface_mesh = surface_mesh_base.copy()
        smp0 = surface_mesh[:,0,:].transpose()
        smp1 = surface_mesh[:,1,:].transpose()
        smp2 = surface_mesh[:,2,:].transpose()
        object_mesh = object_mesh_base.copy()
        # raise ValueError(smp0.shape)

        objp0 = object_mesh[:,0,:].transpose() + np.array([0.001, 0.001, 0.001]).reshape(3,1)
        objp1 = object_mesh[:,1,:].transpose() + np.array([0.001, 0.001, 0.001]).reshape(3,1)
        objp2 = object_mesh[:,2,:].transpose() + np.array([0.001, 0.001, 0.001]).reshape(3,1)           
        result = g.compute(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)
        perim = result[1]
        custom_assert(self, perim, 4*np.sqrt(2)+2*2, base_tol=1e-6)
        helper_test_derivatives_translate_objects_random(self, objp0, objp1, objp2, smp0, smp1, smp2, 10)


class OffsetSphereIntersectedTestCase(unittest.TestCase):
    def test_fixed_intersection(self):
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        test_data_path = self.base_path + r'/inputFiles'
        surface_mesh_base = mesh.Mesh.from_file(test_data_path+'/25mmsphere_reduced.stl').vectors
        object_mesh_base = mesh.Mesh.from_file(test_data_path+'/25mmsphere_reduced.stl').vectors
        smSize = surface_mesh_base[:,0,:].shape[0]
        omSize = object_mesh_base[:,0,:].shape[0]

        surface_mesh = surface_mesh_base.copy()
        smp0 = surface_mesh[:,0,:].transpose()
        smp1 = surface_mesh[:,1,:].transpose()
        smp2 = surface_mesh[:,2,:].transpose()
            
        object_mesh = object_mesh_base.copy()
        offsetdir = np.array([0.,1.,0.])
        sphere_rad = 25 / 2
        offsetmagnitude = 0.8 * sphere_rad
        offsetvec = offsetmagnitude * offsetdir
        offsetvec = offsetvec.reshape((3,1))

        objp0 = object_mesh[:,0,:].transpose() + offsetvec
        objp1 = object_mesh[:,1,:].transpose() + offsetvec
        objp2 = object_mesh[:,2,:].transpose() + offsetvec            
        result = g.compute(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)

        exact_int = 2 * np.pi * sphere_rad * np.sqrt(1 - (offsetmagnitude / 2 / sphere_rad) ** 2)
        custom_assert(self, result[1], 2 * np.pi * sphere_rad * np.sqrt(1 - (offsetmagnitude / 2 / sphere_rad) ** 2), base_tol=5e-3)
        helper_test_derivatives_translate_objects_random(self, objp0, objp1, objp2, smp0, smp1, smp2, 2)
        
    def test_randomly_generated_intersections(self):
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        test_data_path = self.base_path + r'/inputFiles'
        surface_mesh_base = mesh.Mesh.from_file(test_data_path+'/25mmsphere_reduced.stl').vectors
        object_mesh_base = mesh.Mesh.from_file(test_data_path+'/25mmsphere_reduced.stl').vectors
        smSize = surface_mesh_base[:,0,:].shape[0]
        omSize = object_mesh_base[:,0,:].shape[0]

        surface_mesh = surface_mesh_base.copy()
        smp0 = surface_mesh[:,0,:].transpose()
        smp1 = surface_mesh[:,1,:].transpose()
        smp2 = surface_mesh[:,2,:].transpose()

        for i in range(5):
            
            object_mesh = object_mesh_base.copy()
            # generate a random unit vector direction
            offsetdir = np.random.uniform(-1,1,3)
            # offsetdir = np.array([0.,1.,0.])
            offsetdir = offsetdir / np.linalg.norm(offsetdir)
            sphere_rad = 25 / 2
            offsetmagnitude = np.random.uniform(0.1, 3.0 * sphere_rad)
            offsetmagnitude = 0.8 * sphere_rad
            offsetvec = offsetmagnitude * offsetdir
            offsetvec = offsetvec.reshape((3,1))

            objp0 = object_mesh[:,0,:].transpose() + offsetvec
            objp1 = object_mesh[:,1,:].transpose() + offsetvec
            objp2 = object_mesh[:,2,:].transpose() + offsetvec            
            result = g.compute(objp0, objp1, objp2, smp0, smp1, smp2, 0.001, 10)

            if offsetmagnitude < 1.95 * sphere_rad:
                # intersecting
                exact_int = 2 * np.pi * sphere_rad * np.sqrt(1 - (offsetmagnitude / 2 / sphere_rad) ** 2)
                custom_assert(self, result[1], 2 * np.pi * sphere_rad * np.sqrt(1 - (offsetmagnitude / 2 / sphere_rad) ** 2), base_tol=5e-3)
                helper_test_derivatives_translate_objects_random(self, objp0, objp1, objp2, smp0, smp1, smp2,1)

            elif offsetmagnitude > 2.0 * sphere_rad:
                # nonintersecting
                exact_int = 0.0
                custom_assert(self, exact_int, result[1], base_tol=1e-3)

if __name__ == '__main__':
    unittest.main()
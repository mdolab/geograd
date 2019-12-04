import numpy as np 
import unittest
from geograd import geograd as g
from geograd_complex import geograd as gcs
from stl import mesh
import os
h = 1e-15

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

    def test_values(self):
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 10)
        self.assertAlmostEqual(2.33166, result[2], 4)
        self.assertAlmostEqual(0.0, result[1], 10)
        result2 = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300)
        self.assertAlmostEqual(-2.3137489765687533, result2[0], 8)
    
    def test_derivatives_CS(self):
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 10)
        result2 = g.compute_derivs(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 10)
        ks_base = result2[0]
        #print(result[0])
        #print(ks_base)
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
        
        obj_mesh_contrib_indices = np.array([63],dtype=np.int_)
        surf_mesh_contrib_indices = np.array([1],dtype=np.int_)

        A1_pts = np.take(self.objp0, obj_mesh_contrib_indices, axis=1)
        B1_pts = np.take(self.objp1, obj_mesh_contrib_indices, axis=1)
        C1_pts = np.take(self.objp2, obj_mesh_contrib_indices, axis=1)
        p1_to_perturb = np.unique(np.hstack([A1_pts, B1_pts, C1_pts]), axis=1)

        tot_counter = 0
        nonzero_grad_md = False
        nonzero_grad_ks = False

        #first perturb surface mesh points
        for i in range(p1_to_perturb.shape[1]):
            # search each input for a particular point and get the indices where that point appears
            vec_to_search_for = p1_to_perturb[:,i].reshape((3,1))
            # print('Search vec: '+str(vec_to_search_for))
            A1_indices = np.argwhere(np.all((self.objp0-vec_to_search_for)==0, axis=0))
            B1_indices = np.argwhere(np.all((self.objp1-vec_to_search_for)==0, axis=0))
            C1_indices = np.argwhere(np.all((self.objp2-vec_to_search_for)==0, axis=0))

        #     # 4) do a complex step increment to each x, y, and z component for each point at all locations in p and q where it appears
        #     # this alters the mesh in a self-consistent fashion same as pyGeo would
            for k in range(3):
                exact_grad_ks = 0
                A1_alt = self.objp0.copy().astype(np.complex64)
                for j in A1_indices:
                    jind = j[0]
                    A1_alt[k,jind] = A1_alt[k,jind] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += A1_grad[k,jind]
                B1_alt = self.objp1.copy().astype(np.complex64)
                for j in B1_indices:
                    jind = j[0]
                    B1_alt[k,jind] = B1_alt[k,jind] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += B1_grad[k,jind]
                C1_alt = self.objp2.copy().astype(np.complex64)
                for j in C1_indices:
                    jind = j[0]
                    C1_alt[k,jind] = C1_alt[k,jind] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += C1_grad[k,jind]

                resultcs = gcs.compute(A1_alt, B1_alt, C1_alt, self.smp0, self.smp1, self.smp2, 1.0, 10)
                resultcs2 = gcs.compute(A1_alt, B1_alt, C1_alt, self.smp0, self.smp1, self.smp2, resultcs[2], 10)
                ks_cs = resultcs2[0]
                gradcs_ks = np.imag((ks_cs - ks_base)) / cseps

                print('Exact: '+str(exact_grad_ks))
                print('cs: '+str(gradcs_ks))
        #         # warnings.warn('Surf mesh CS: '+str(gradcs_ks)+' Exact: '+str(exact_grad_ks))
        #         tot_counter += 1
        #         if np.abs(gradcs_ks) > 1e-3:
        #             nonzero_grad_ks = True

        # self.assertTrue(tot_counter > 0, msg='If tot_counter remains zero there is some issue finding close points using the numpy slicing and searching')
        # self.assertTrue(nonzero_grad_ks, msg='Check to make sure at least one gradient checked is actually nonzero')


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

    def test_values(self):
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 10)
        self.assertAlmostEqual(0.0178387, result[2], 6)
        self.assertAlmostEqual(0.0, result[1], 10)
        result2 = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300)
        self.assertAlmostEqual(-0.013197699692565058, result2[0], 8)


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

    def test_values(self):
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 10)
        self.assertAlmostEqual(0.02212236, result[2], 6)
        self.assertAlmostEqual(0.0, result[1], 10)
        result2 = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300)
        self.assertAlmostEqual(-0.015186799790516424, result2[0], 8)
        result3 = g.compute_derivs(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300)
        # md indices 104, 52
        # get the ABC1 points at index 104 and maybe a couple of random others
        # for each of the points:
        # find all the points in A1, B1, C1 that match

        # get all the points at 52
    
    def test_derivatives_CS(self):
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 300)
        result2 = g.compute_derivs(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300)
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
        
        obj_mesh_contrib_indices = np.array([104],dtype=np.int_)
        surf_mesh_contrib_indices = np.array([52],dtype=np.int_)

        A1_pts = np.take(self.objp0, obj_mesh_contrib_indices, axis=1)
        B1_pts = np.take(self.objp1, obj_mesh_contrib_indices, axis=1)
        C1_pts = np.take(self.objp2, obj_mesh_contrib_indices, axis=1)
        p1_to_perturb = np.unique(np.hstack([A1_pts, B1_pts, C1_pts]), axis=1)

        tot_counter = 0
        nonzero_grad_md = False
        nonzero_grad_ks = False

        #first perturb surface mesh points
        for i in range(p1_to_perturb.shape[1]):
            # search each input for a particular point and get the indices where that point appears
            vec_to_search_for = p1_to_perturb[:,i].reshape((3,1))
            # print('Search vec: '+str(vec_to_search_for))
            A1_indices = np.argwhere(np.all((self.objp0-vec_to_search_for)==0, axis=0))
            B1_indices = np.argwhere(np.all((self.objp1-vec_to_search_for)==0, axis=0))
            C1_indices = np.argwhere(np.all((self.objp2-vec_to_search_for)==0, axis=0))

        #     # 4) do a complex step increment to each x, y, and z component for each point at all locations in p and q where it appears
        #     # this alters the mesh in a self-consistent fashion same as pyGeo would
            if True:
                for k in range(3):
                    exact_grad_ks = 0
                    A1_alt = self.objp0.copy().astype(np.complex64)
                    for j in A1_indices:
                        jind = j[0]
                        A1_alt[k,jind] = A1_alt[k,jind] + cseps*1.0j
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += A1_grad[k,jind]
                    B1_alt = self.objp1.copy().astype(np.complex64)
                    for j in B1_indices:
                        jind = j[0]
                        B1_alt[k,jind] = B1_alt[k,jind] + cseps*1.0j
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += B1_grad[k,jind]
                    C1_alt = self.objp2.copy().astype(np.complex64)
                    for j in C1_indices:
                        jind = j[0]
                        C1_alt[k,jind] = C1_alt[k,jind] + cseps*1.0j
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += C1_grad[k,jind]

                    resultcs = gcs.compute(A1_alt, B1_alt, C1_alt, self.smp0, self.smp1, self.smp2, 1.0, 300)
                    resultcs2 = gcs.compute(A1_alt, B1_alt, C1_alt, self.smp0, self.smp1, self.smp2, resultcs[2], 300)
                    ks_cs = resultcs2[0]
                    gradcs_ks = np.imag((ks_cs - ks_base)) / cseps

                    print('Exact: '+str(exact_grad_ks))
                    print('cs: '+str(gradcs_ks))
            else:
                for k in range(3):
                    exact_grad_ks = 0
                    A1_alt = self.objp0.copy()
                    for j in A1_indices:
                        jind = j[0]
                        A1_alt[k,jind] = A1_alt[k,jind] + fdeps
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += A1_grad[k,jind]
                    B1_alt = self.objp1.copy()
                    for j in B1_indices:
                        jind = j[0]
                        B1_alt[k,jind] = B1_alt[k,jind] + fdeps
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += B1_grad[k,jind]
                    C1_alt = self.objp2.copy()
                    for j in C1_indices:
                        jind = j[0]
                        C1_alt[k,jind] = C1_alt[k,jind] + fdeps
                        # 5) sum the contributions from each point from the analytic gradient and compare
                        exact_grad_ks += C1_grad[k,jind]

                    resultfd = g.compute(A1_alt, B1_alt, C1_alt, self.smp0, self.smp1, self.smp2, 1.0, 300)
                    resultfd2 = g.compute(A1_alt, B1_alt, C1_alt, self.smp0, self.smp1, self.smp2, resultfd[2], 300)
                    ks_fd = resultfd2[0]
                    gradfd_ks = (ks_fd - ks_base) / fdeps

                    print('Exact: '+str(exact_grad_ks))
                    print('fd: '+str(gradfd_ks))
        #         # warnings.warn('Surf mesh CS: '+str(gradcs_ks)+' Exact: '+str(exact_grad_ks))
        #         tot_counter += 1
        #         if np.abs(gradcs_ks) > 1e-3:
        #             nonzero_grad_ks = True

        # self.assertTrue(tot_counter > 0, msg='If tot_counter remains zero there is some issue finding close points using the numpy slicing and searching')
        # self.assertTrue(nonzero_grad_ks, msg='Check to make sure at least one gradient checked is actually nonzero')

# test the other two blobs

if __name__ == '__main__':
    unittest.main()
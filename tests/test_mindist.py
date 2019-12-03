import geograd
import tensorflow as tf
tf.enable_eager_execution()
from geograd.minimum_distance import mindist
import pandas as pd
import numpy as np
from openmdao.utils.assert_utils import assert_rel_error
from numpy.testing import assert_almost_equal
import unittest
import os
from stl import mesh
import warnings

def custom_assert(self, truth, approx, base_tol=1e-6):
    if np.abs(truth) > 0.1:
        assert_rel_error(self, truth, approx, tolerance=base_tol)
    elif np.abs(truth) > 1e-4:
        assert_rel_error(self, truth, approx, tolerance=base_tol*10)
    elif np.abs(truth) > 5e-7:
        assert_rel_error(self, truth, approx, tolerance=base_tol*100)
    else:
        assert_almost_equal(truth, approx, decimal=7)

class MinDistSTLTestCase1(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.getcwd()+'/geograd/tests'):
            test_data_path = os.getcwd()+'/geograd/tests/data'
        elif os.path.isdir(os.getcwd()+'/tests'):
            test_data_path = os.getcwd()+'/tests/data'
        else:
            test_data_path = os.getcwd()+'/data'

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

    def test_values_KS(self):
        KS, min_dist, blah1, blah2 = mindist(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=False, constraint_type='KS')
        assert_rel_error(self, 2.33166, min_dist, tolerance=1e-4)
        assert_rel_error(self, -2.31435, KS, tolerance=1e-4)

    def test_identical_KS(self):
        KS1, min_dist1, blah1, blah2 = mindist(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=False, constraint_type='KS')
        KS2, min_dist2, blah1, blah2 = mindist(self.smp0.reshape(self.smSize, 1, 3), self.smp1.reshape(self.smSize, 1, 3), self.smp2.reshape(self.smSize, 1, 3), 
                                               self.objp0.reshape(1, self.omSize, 3), self.objp1.reshape(1, self.omSize, 3), self.objp2.reshape(1, self.omSize, 3), 
                                               dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=False, constraint_type='KS')

        assert_rel_error(self, 2.33166, min_dist1, tolerance=1e-4)
        assert_rel_error(self, -2.31435, KS1, tolerance=1e-4)
        assert_almost_equal(min_dist1-min_dist2, 0., decimal=15)
        assert_almost_equal(KS1-KS2, 0., decimal=9)


    def test_values_DIE(self):
        DIE, min_dist, blah1, blah2 = mindist(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 20, batch_size=20, complexify=False, constraint_type='DIE')
        assert_rel_error(self, 2.33166, min_dist, tolerance=1e-4)
        assert_rel_error(self, -2.41516, DIE, tolerance=1e-4)

class MinDistSTLTestCase2(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.getcwd()+'/geograd/tests'):
            test_data_path = os.getcwd()+'/geograd/tests/data'
        elif os.path.isdir(os.getcwd()+'/tests'):
            test_data_path = os.getcwd()+'/tests/data'
        else:
            test_data_path = os.getcwd()+'/data'

        # the 'surface'  mesh is a the same blob file offset by 3 units in the y direction
        surface_mesh = mesh.Mesh.from_file(test_data_path+'/blob_offset_y3_line_line.stl').vectors
        object_mesh = mesh.Mesh.from_file(test_data_path+'/blob_line_line.stl').vectors
        self.smSize = surface_mesh[:,0,:].shape[0]
        self.smp0 = surface_mesh[:,0,:].reshape(1, self.smSize, 3)
        self.smp1 = surface_mesh[:,1,:].reshape(1, self.smSize, 3)
        self.smp2 = surface_mesh[:,2,:].reshape(1, self.smSize, 3)

        self.omSize = object_mesh[:,0,:].shape[0]
        self.objp0 = object_mesh[:,0,:].reshape(self.omSize, 1, 3)
        self.objp1 = object_mesh[:,1,:].reshape(self.omSize, 1, 3)
        self.objp2 = object_mesh[:,2,:].reshape(self.omSize, 1, 3)

    def test_values_KS(self):
        KS, min_dist, blah1, blah2 = mindist(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=False, constraint_type='KS')
        assert_rel_error(self, 0.0178387, min_dist, tolerance=1e-4)
        assert_rel_error(self, -0.01319848588256959, KS, tolerance=1e-4)

    def test_values_DIE(self):
        DIE, min_dist, blah1, blah2 = mindist(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 20, batch_size=20, complexify=False, constraint_type='DIE')
        assert_rel_error(self, 0.0178387, min_dist, tolerance=1e-4)
        assert_rel_error(self, -0.12958106798770888, DIE, tolerance=1e-4)

    def test_identical_KS(self):
        KS1, min_dist1, blah1, blah2 = mindist(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=False, constraint_type='KS')
        KS2, min_dist2, blah1, blah2 = mindist(self.smp0.reshape(self.smSize, 1, 3), self.smp1.reshape(self.smSize, 1, 3), self.smp2.reshape(self.smSize, 1, 3), 
                                               self.objp0.reshape(1, self.omSize, 3), self.objp1.reshape(1, self.omSize, 3), self.objp2.reshape(1, self.omSize, 3), 
                                               dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=False, constraint_type='KS')

        assert_rel_error(self, 0.0178387, min_dist1, tolerance=1e-4)
        assert_rel_error(self, -0.01319848588256959, KS1, tolerance=1e-4)
        assert_almost_equal(min_dist1-min_dist2, 0., decimal=15)
        assert_almost_equal(KS1-KS2, 0., decimal=5)

    def test_CS_derivs_KS(self):
        ks, min_dist, all_grads, blah2, pairwise_dists, PT_closer, stacks= mindist(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=False, return_extras=True)

        A1_grad = all_grads[0]
        B1_grad = all_grads[1]
        C1_grad = all_grads[2]
        A2_grad = all_grads[3]
        B2_grad = all_grads[4]
        C2_grad = all_grads[5]

        # because of the way STL files are stored, mesh vertices appear multiple times in the ABC matrices.
        # they need to be manipulated at the same time or the minimum distance gradients won't be correct

        cseps = 1.0e-10
        fdeps = 5.0e-3
        #raise ValueError(np.argwhere(np.abs(A1_grad)>1e-6))
        # find the indices of point which alter the minimum distance and manipulate them all together
        equal_min_indices= np.argwhere(pairwise_dists==min_dist)
        # also search a random couple of points not on the min dist path but with finite grads
        obj_mesh_contrib_indices = np.append(np.unique(equal_min_indices[:,0]),[63+A1_grad.shape[0]])
        surf_mesh_contrib_indices = np.append(np.unique(equal_min_indices[:,1]),[50])

        if PT_closer:
            # the point-triangle tests are governing the minimum distance computation
            obj_md_pts = np.unique(np.take(stacks[0], obj_mesh_contrib_indices, axis=0),axis=0)
            mesh_A_pts = np.take(self.smp0, surf_mesh_contrib_indices, axis=1)
            mesh_B_pts = np.take(self.smp1, surf_mesh_contrib_indices, axis=1)
            mesh_C_pts = np.take(self.smp2, surf_mesh_contrib_indices, axis=1)
            mesh_md_pts = np.unique(np.hstack([mesh_A_pts, mesh_B_pts, mesh_C_pts]),axis=1)

        else:
            # the line-line tests are governing the minimum distance computation
            obj_p_pts = np.take(stacks[0], obj_mesh_contrib_indices, axis=0)
            obj_q_pts = np.take(stacks[1], obj_mesh_contrib_indices, axis=0)
            obj_md_pts = np.unique(np.vstack([obj_p_pts, obj_q_pts]),axis=0)
            mesh_p_pts = np.take(stacks[2], surf_mesh_contrib_indices, axis=1)
            mesh_q_pts = np.take(stacks[3], surf_mesh_contrib_indices, axis=1)
            mesh_md_pts = np.unique(np.hstack([mesh_p_pts, mesh_q_pts]),axis=1)

        tot_counter = 0
        nonzero_grad_md = False
        nonzero_grad_ks = False

        #first perturb surface mesh points
        for i in range(mesh_md_pts.shape[1]):
            # search each input for a particular point and get the indices where that point appears
            vec_to_search_for = np.squeeze(mesh_md_pts[:,i,:])
            A2_indices = np.argwhere(np.all((self.smp0-vec_to_search_for)==0, axis=2))[:,1]
            B2_indices = np.argwhere(np.all((self.smp1-vec_to_search_for)==0, axis=2))[:,1]
            C2_indices = np.argwhere(np.all((self.smp2-vec_to_search_for)==0, axis=2))[:,1]

            # 4) do a complex step increment to each x, y, and z component for each point at all locations in p and q where it appears
            # this alters the mesh in a self-consistent fashion same as pyGeo would
            for k in range(3):
                exact_grad_ks = 0

                mesh_A_alt = self.smp0.copy().astype(np.complex64)
                for j in A2_indices:
                    mesh_A_alt[0,j,k] = mesh_A_alt[0,j,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += A2_grad[j,k]

                mesh_B_alt = self.smp1.copy().astype(np.complex64)
                for j in B2_indices:
                    mesh_B_alt[0,j,k] = mesh_B_alt[0,j,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += B2_grad[j,k]

                mesh_C_alt = self.smp2.copy().astype(np.complex64)
                for j in C2_indices:
                    mesh_C_alt[0,j,k] = mesh_C_alt[0,j,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += C2_grad[j,k]

                kscs, blah0 = mindist(self.objp0, self.objp1, self.objp2, mesh_A_alt, mesh_B_alt, mesh_C_alt, dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=True)
                gradcs_ks = np.imag((kscs - ks)) / cseps

                custom_assert(self, exact_grad_ks, gradcs_ks, base_tol=1e-6)
                # warnings.warn('Surf mesh CS: '+str(gradcs_ks)+' Exact: '+str(exact_grad_ks))
                tot_counter += 1
                if np.abs(gradcs_ks) > 1e-3:
                    nonzero_grad_ks = True

        self.assertTrue(tot_counter > 0, msg='If tot_counter remains zero there is some issue finding close points using the numpy slicing and searching')
        self.assertTrue(nonzero_grad_ks, msg='Check to make sure at least one gradient checked is actually nonzero')

        # 3) locate the unique xyz points in the p and q matrix on each side
        #in general the point of interest will occur more times than in step 2 since it is connected to non critical lines
        # for an object mesh test:
        tot_counter = 0
        nonzero_grad_ks = False

        for i in range(obj_md_pts.shape[0]):
            vec_to_search_for = np.squeeze(obj_md_pts[i,:,:])
            A1_indices = np.argwhere(np.all((self.objp0-vec_to_search_for)==0, axis=2))[:,0]
            B1_indices = np.argwhere(np.all((self.objp1-vec_to_search_for)==0, axis=2))[:,0]
            C1_indices = np.argwhere(np.all((self.objp2-vec_to_search_for)==0, axis=2))[:,0]

            # 4) do a complex step increment to each x, y, and z component for each point at all locations in p and q where it appears
            # this alters the mesh in a self-consistent fashion same as pyGeo would
            for k in range(3):
                exact_grad_ks = 0
                exact_grad_md = 0

                obj_A_alt = self.objp0.copy().astype(np.complex64)
                for j in A1_indices:
                    obj_A_alt[j,0,k] = obj_A_alt[j,0,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += A1_grad[j,k]

                obj_B_alt = self.objp1.copy().astype(np.complex64)
                for j in B1_indices:
                    obj_B_alt[j,0,k] = obj_B_alt[j,0,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += B1_grad[j,k]

                obj_C_alt = self.objp2.copy().astype(np.complex64)
                for j in C1_indices:
                    obj_C_alt[j,0,k] = obj_C_alt[j,0,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += C1_grad[j,k]

                kscs, blah0 = mindist(obj_A_alt, obj_B_alt, obj_C_alt, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=True)
                gradcs_ks = np.imag((kscs - ks)) / cseps
                custom_assert(self, exact_grad_ks, gradcs_ks, base_tol=1e-6)
                # warnings.warn('Obj mesh CS: '+str(gradcs_ks)+' Exact: '+str(exact_grad_ks) + ' Vec: '+str(vec_to_search_for)+' k: '+str(k),UserWarning)

                tot_counter += 1
                if np.abs(gradcs_ks) > 1e-3:
                    nonzero_grad_ks = True

        self.assertTrue(tot_counter > 0, msg='If tot_counter remains zero there is some issue finding close points using the numpy slicing and searching')
        self.assertTrue(nonzero_grad_ks, msg='Check to make sure at least one gradient checked is actually nonzero')


class MinDistSTLTestCase3(unittest.TestCase):
    def setUp(self):
        if os.path.isdir(os.getcwd()+'/geograd/tests'):
            test_data_path = os.getcwd()+'/geograd/tests/data'
        elif os.path.isdir(os.getcwd()+'/tests'):
            test_data_path = os.getcwd()+'/tests/data'
        else:
            test_data_path = os.getcwd()+'/data'

        # the 'surface'  mesh is a the same blob file offset by 3 units in the y direction
        surface_mesh = mesh.Mesh.from_file(test_data_path+'/blob_offset_y3_point_tri.stl').vectors
        object_mesh = mesh.Mesh.from_file(test_data_path+'/blob_point_tri.stl').vectors
        self.smSize = surface_mesh[:,0,:].shape[0]
        self.smp0 = surface_mesh[:,0,:].reshape(1, self.smSize, 3)
        self.smp1 = surface_mesh[:,1,:].reshape(1, self.smSize, 3)
        self.smp2 = surface_mesh[:,2,:].reshape(1, self.smSize, 3)

        self.omSize = object_mesh[:,0,:].shape[0]
        self.objp0 = object_mesh[:,0,:].reshape(self.omSize, 1, 3)
        self.objp1 = object_mesh[:,1,:].reshape(self.omSize, 1, 3)
        self.objp2 = object_mesh[:,2,:].reshape(self.omSize, 1, 3)

    def test_values_KS(self):
        KS, min_dist, blah1, blah2 = mindist(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=False, constraint_type='KS')
        assert_rel_error(self, 0.02212236, min_dist, tolerance=1e-4)
        assert_rel_error(self, -0.015186954228556992, KS, tolerance=1e-4)

    def test_values_DIE(self):
        DIE, min_dist, blah1, blah2 = mindist(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 20, batch_size=20, complexify=False, constraint_type='DIE')
        assert_rel_error(self, 0.02212236, min_dist, tolerance=1e-4)
        assert_rel_error(self, -0.11771683644963829, DIE, tolerance=1e-4)

    def test_CS_derivs_KS(self):
        ks, min_dist, all_grads, blah2, pairwise_dists, PT_closer, stacks= mindist(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=False, return_extras=True)

        A1_grad = all_grads[0]
        B1_grad = all_grads[1]
        C1_grad = all_grads[2]
        A2_grad = all_grads[3]
        B2_grad = all_grads[4]
        C2_grad = all_grads[5]

        # because of the way STL files are stored, mesh vertices appear multiple times in the ABC matrices.
        # they need to be manipulated at the same time or the minimum distance gradients won't be correct

        cseps = 1.0e-10
        fdeps = 5.0e-3
        #raise ValueError(np.argwhere(np.abs(A1_grad)>1e-6))
        # find the indices of point which alter the minimum distance and manipulate them all together
        equal_min_indices= np.argwhere(pairwise_dists==min_dist)
        # also search a random couple of points not on the min dist path but with finite grads
        obj_mesh_contrib_indices = np.unique(equal_min_indices[:,0])
        surf_mesh_contrib_indices = np.unique(equal_min_indices[:,1])

        if PT_closer:
            # the point-triangle tests are governing the minimum distance computation
            obj_md_pts = np.unique(np.take(stacks[0], obj_mesh_contrib_indices, axis=0),axis=0)
            mesh_A_pts = np.take(self.smp0, surf_mesh_contrib_indices, axis=1)
            mesh_B_pts = np.take(self.smp1, surf_mesh_contrib_indices, axis=1)
            mesh_C_pts = np.take(self.smp2, surf_mesh_contrib_indices, axis=1)
            mesh_md_pts = np.unique(np.hstack([mesh_A_pts, mesh_B_pts, mesh_C_pts]),axis=1)

        else:
            # the line-line tests are governing the minimum distance computation
            obj_p_pts = np.take(stacks[0], obj_mesh_contrib_indices, axis=0)
            obj_q_pts = np.take(stacks[1], obj_mesh_contrib_indices, axis=0)
            obj_md_pts = np.unique(np.vstack([obj_p_pts, obj_q_pts]),axis=0)
            mesh_p_pts = np.take(stacks[2], surf_mesh_contrib_indices, axis=1)
            mesh_q_pts = np.take(stacks[3], surf_mesh_contrib_indices, axis=1)
            mesh_md_pts = np.unique(np.hstack([mesh_p_pts, mesh_q_pts]),axis=1)

        tot_counter = 0
        nonzero_grad_md = False
        nonzero_grad_ks = False

        #first perturb surface mesh points
        for i in range(mesh_md_pts.shape[1]):
            # search each input for a particular point and get the indices where that point appears
            vec_to_search_for = np.squeeze(mesh_md_pts[:,i,:])
            A2_indices = np.argwhere(np.all((self.smp0-vec_to_search_for)==0, axis=2))[:,1]
            B2_indices = np.argwhere(np.all((self.smp1-vec_to_search_for)==0, axis=2))[:,1]
            C2_indices = np.argwhere(np.all((self.smp2-vec_to_search_for)==0, axis=2))[:,1]

            # 4) do a complex step increment to each x, y, and z component for each point at all locations in p and q where it appears
            # this alters the mesh in a self-consistent fashion same as pyGeo would
            for k in range(3):
                exact_grad_ks = 0

                mesh_A_alt = self.smp0.copy().astype(np.complex64)
                for j in A2_indices:
                    mesh_A_alt[0,j,k] = mesh_A_alt[0,j,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += A2_grad[j,k]

                mesh_B_alt = self.smp1.copy().astype(np.complex64)
                for j in B2_indices:
                    mesh_B_alt[0,j,k] = mesh_B_alt[0,j,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += B2_grad[j,k]

                mesh_C_alt = self.smp2.copy().astype(np.complex64)
                for j in C2_indices:
                    mesh_C_alt[0,j,k] = mesh_C_alt[0,j,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += C2_grad[j,k]

                kscs, blah0 = mindist(self.objp0, self.objp1, self.objp2, mesh_A_alt, mesh_B_alt, mesh_C_alt, dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=True)
                gradcs_ks = np.imag((kscs - ks)) / cseps

                custom_assert(self, exact_grad_ks, gradcs_ks, base_tol=1e-6)
                # warnings.warn('Surf mesh CS: '+str(gradcs_ks)+' Exact: '+str(exact_grad_ks))
                tot_counter += 1
                if np.abs(gradcs_ks) > 1e-3:
                    nonzero_grad_ks = True

        self.assertTrue(tot_counter > 0, msg='If tot_counter remains zero there is some issue finding close points using the numpy slicing and searching')
        self.assertTrue(nonzero_grad_ks, msg='Check to make sure at least one gradient checked is actually nonzero')

        # 3) locate the unique xyz points in the p and q matrix on each side
        #in general the point of interest will occur more times than in step 2 since it is connected to non critical lines
        # for an object mesh test:
        tot_counter = 0
        nonzero_grad_ks = False

        for i in range(obj_md_pts.shape[0]):
            vec_to_search_for = np.squeeze(obj_md_pts[i,:,:])
            A1_indices = np.argwhere(np.all((self.objp0-vec_to_search_for)==0, axis=2))[:,0]
            B1_indices = np.argwhere(np.all((self.objp1-vec_to_search_for)==0, axis=2))[:,0]
            C1_indices = np.argwhere(np.all((self.objp2-vec_to_search_for)==0, axis=2))[:,0]

            # 4) do a complex step increment to each x, y, and z component for each point at all locations in p and q where it appears
            # this alters the mesh in a self-consistent fashion same as pyGeo would
            for k in range(3):
                exact_grad_ks = 0
                exact_grad_md = 0

                obj_A_alt = self.objp0.copy().astype(np.complex64)
                for j in A1_indices:
                    obj_A_alt[j,0,k] = obj_A_alt[j,0,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += A1_grad[j,k]

                obj_B_alt = self.objp1.copy().astype(np.complex64)
                for j in B1_indices:
                    obj_B_alt[j,0,k] = obj_B_alt[j,0,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += B1_grad[j,k]

                obj_C_alt = self.objp2.copy().astype(np.complex64)
                for j in C1_indices:
                    obj_C_alt[j,0,k] = obj_C_alt[j,0,k] + cseps*1.0j
                    # 5) sum the contributions from each point from the analytic gradient and compare
                    exact_grad_ks += C1_grad[j,k]

                kscs, blah0 = mindist(obj_A_alt, obj_B_alt, obj_C_alt, self.smp0, self.smp1, self.smp2, dist_tol=0.0, adapt_rho=False, rho_c = 300, batch_size=20, complexify=True)
                gradcs_ks = np.imag((kscs - ks)) / cseps
                custom_assert(self, exact_grad_ks, gradcs_ks, base_tol=1e-6)
                # warnings.warn('Obj mesh CS: '+str(gradcs_ks)+' Exact: '+str(exact_grad_ks) + ' Vec: '+str(vec_to_search_for)+' k: '+str(k),UserWarning)

                tot_counter += 1
                if np.abs(gradcs_ks) > 1e-3:
                    nonzero_grad_ks = True

        self.assertTrue(tot_counter > 0, msg='If tot_counter remains zero there is some issue finding close points using the numpy slicing and searching')
        self.assertTrue(nonzero_grad_ks, msg='Check to make sure at least one gradient checked is actually nonzero')

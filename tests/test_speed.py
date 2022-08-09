import numpy as np 
import unittest
from geograd import geograd_parallel as g
from geograd import geograd_parallel_complex as gcs
from stl import mesh
import os
h = 1e-15
from mpi4py import MPI
import time

def custom_assert(self, truth, approx, base_tol=1e-7):
    if np.abs(truth) > 0.1:
        assert_rel_error(self, truth, approx, tolerance=base_tol)
    elif np.abs(truth) > 1e-4:
        assert_rel_error(self, truth, approx, tolerance=base_tol*10)
    elif np.abs(truth) > 5e-9:
        assert_rel_error(self, truth, approx, tolerance=base_tol*100)
    else:
        assert_almost_equal(truth, approx, decimal=7)

# class MinDistSTLTestCase1(unittest.TestCase):
#     def setUp(self):
#         if os.path.isdir(os.getcwd()+'/tests'):
#             test_data_path = os.getcwd()+'/tests'
#         else:
#             raise IOError

#         # the 'surface'  mesh is a the same blob file offset by 3 units in the y direction
#         surface_mesh = mesh.Mesh.from_file(test_data_path+'/blob_offset_y3.stl').vectors
#         object_mesh = mesh.Mesh.from_file(test_data_path+'/blob.stl').vectors
#         self.smSize = surface_mesh[:,0,:].shape[0]
#         self.smp0 = surface_mesh[:,0,:].transpose()
#         self.smp1 = surface_mesh[:,1,:].transpose()
#         self.smp2 = surface_mesh[:,2,:].transpose()

#         self.omSize = object_mesh[:,0,:].shape[0]
#         self.objp0 = object_mesh[:,0,:].transpose()
#         self.objp1 = object_mesh[:,1,:].transpose()
#         self.objp2 = object_mesh[:,2,:].transpose()

#     def test_values(self):
#         start=time.time()
#         for i in range(100):
#             result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0, 10, 2.0)
#         end = time.time()
#         if MPI.COMM_WORLD.rank == 0:
#             print('Elapsed time for 100 runs of analysis only: '+str(end-start))

#         self.assertAlmostEqual(2.33166, result[2], 4)
#         self.assertAlmostEqual(0.0, result[1], 10)


#         start = time.time()
#         for i in range(100):
#             result2 = g.compute_derivs(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, result[2], 300, 2.0)
#         end = time.time()
#         if MPI.COMM_WORLD.rank == 0:
#             print('Elapsed time for 100 runs with derivatives: '+str(end-start))


class BWBTestCase(unittest.TestCase):
    def setUp(self):
        self.base_path = os.path.dirname(os.path.abspath(__file__))
        test_data_path = self.base_path + r'/inputFiles'

        # the 'surface'  mesh is a the same blob file offset by 3 units in the y direction
        surface_mesh = mesh.Mesh.from_file(test_data_path+'/bwb.stl').vectors
        object_mesh = mesh.Mesh.from_file(test_data_path+'/blob_bwb_wing.stl').vectors
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
        transl_vec = np.linspace(0,28.0, 50)
        result = g.compute(self.smp0, self.smp1, self.smp2, self.objp0, self.objp1, self.objp2, 1.0, 10, self.maxdim, MPI.COMM_WORLD.py2f())
        self.assertAlmostEqual(0.1698771104979689, result[2], 4)
        self.assertAlmostEqual(0.0, result[1], 10)
        result = g.compute(self.smp0, self.smp1, self.smp2, self.objp0, self.objp1, self.objp2, result[2], 10, self.maxdim, MPI.COMM_WORLD.py2f())
        self.assertAlmostEqual(0.3466089344854115, result[0], 4)
        loop_time = 0
        unbalance = 0
        start=time.time()
        for i in range(50):
            offset = np.array([0.0,  0.0, -transl_vec[i]]).reshape(3,1)
            result = g.compute(self.smp0, self.smp1, self.smp2, self.objp0+offset, self.objp1+offset, self.objp2+offset, 1.0, 10, self.maxdim, MPI.COMM_WORLD.py2f())
            loop_time = loop_time+result[3]
            unbalance += result[4]
        end = time.time()
        if MPI.COMM_WORLD.rank == 0:
            print('Elapsed time for 50 runs of analysis only: '+str(end-start))
            print('Loop time: '+str(loop_time))
            print('Avg loop utilization: '+str(unbalance/50)+' percent')
        loop_time = 0
        unbalance = 0
        start = time.time()
        for i in range(50):
            offset = np.array([0.0,  0.0, -transl_vec[i]]).reshape(3,1)
            # result = g.compute(self.smp0, self.smp1, self.smp2, self.objp0+offset, self.objp1+offset, self.objp2+offset, 1.0, 10, self.maxdim)
            result = g.compute_derivs(self.smp0, self.smp1, self.smp2, self.objp0+offset, self.objp1+offset, self.objp2+offset, 0.10, 10, self.maxdim, MPI.COMM_WORLD.py2f())
            loop_time = loop_time+result[3]
            unbalance += result[4]
        end = time.time()
        if MPI.COMM_WORLD.rank == 0:
            print('Elapsed time for 50 runs with derivatives: '+str(end-start))
            print('Loop time: '+str(loop_time))
            print('Avg loop utilization: '+str(unbalance/50)+' percent')


if __name__ == '__main__':
    unittest.main()
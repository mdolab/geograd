import numpy as np 
import unittest
from geograd import geograd as g
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
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0)
        self.assertAlmostEqual(result[1], 0.5*4, 10)

    def test_intersect_2_in_1(self):
        a2 = np.repeat(np.array([0.75, 0.5, -2.0]),2).reshape(3,2)
        b2 = np.repeat(np.array([-0.25, 0.5, -2.0]),2).reshape(3,2)
        c2 = np.repeat(np.array([0.25, 0.5, 1.0]),2).reshape(3,2)
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0)
        self.assertAlmostEqual(result[1], 1/3*4, 10)

    def test_intersect_right_edge(self):
        a2 = np.repeat(np.array([1.0, 0.5, -2.0]),2).reshape(3,2)
        b2 = np.repeat(np.array([0.0, 0.5, -2.0]),2).reshape(3,2)
        c2 = np.repeat(np.array([0.5, 0.5, 1.0]),2).reshape(3,2)
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0)
        self.assertAlmostEqual(result[1], 1/6*4, 10)

    def test_intersect_left_edge(self):
        a2 = np.repeat(np.array([0.5, 0.5, -2.0]),2).reshape(3,2)
        b2 = np.repeat(np.array([-0.5, 0.5, -2.0]),2).reshape(3,2)
        c2 = np.repeat(np.array([0.0, 0.5, 1.0]),2).reshape(3,2)
        result = g.compute(self.a1, self.b1, self.c1, a2, b2, c2, 1.0)
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
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0)
        self.assertAlmostEqual(2.33166**2, result[2], 4)
        self.assertAlmostEqual(0.0, result[1], 10)

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
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0)
        self.assertAlmostEqual(0.0178387**2, result[2], 6)
        self.assertAlmostEqual(0.0, result[1], 10)



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
        result = g.compute(self.objp0, self.objp1, self.objp2, self.smp0, self.smp1, self.smp2, 1.0)
        self.assertAlmostEqual(0.02212236**2, result[2], 6)
        self.assertAlmostEqual(0.0, result[1], 10)


# test the other two blobs

if __name__ == '__main__':
    unittest.main()
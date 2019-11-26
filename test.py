import numpy as np 
import unittest
from triangles import triangles as t 

class TestPointTri(unittest.TestCase):
    def setUp(self):
        self.a = np.array([0.0, 0.0, 0.0])
        self.b = np.array([1.0, 0.0, 0.0])
        self.c = np.array([0.0, 1.0, 0.0])
    
    def test_interior_plane(self):
        p = np.array([0.25, 0.25, 0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.0, 10)

    def test_interior_nonplane(self):
        p = np.array([0.25, 0.25, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.25, 10)

    def test_edges_plane(self):
        # edge ab
        p = np.array([0.5, -0.5, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.25, 10)
        # edge bc
        p = np.array([1.0, 1.0, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)
        # edge ac
        p = np.array([-0.5, 0.5, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.25, 10)

    def test_edges_nonplane(self):
        # edge ab
        p = np.array([0.5, -0.5, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)
        # edge bc
        p = np.array([1.0, 1.0, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.75, 10)
        # edge ac
        p = np.array([-0.5, 0.5, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)

    def test_vertices_plane(self):
        # vertex a
        p = np.array([-0.5, -0.5, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)
        # vertex b
        p = np.array([1.5, 0.0, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.25, 10)
        # vertex c
        p = np.array([0.0, 1.5, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.25, 10)

    def test_vertices_nonplane(self):
        # vertex a
        p = np.array([-0.5, -0.5, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.75, 10)
        # vertex b
        p = np.array([1.5, 0.0, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)
        # vertex c
        p = np.array([0.0, 1.5, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)

class TestLineLine(unittest.TestCase):
    def setUp(self):
        self.p1 = np.array([0.0,0.0,0.0])
        self.q1 = np.array([1.0, 1.0, 0.0])

    def test_intersect(self):
        # clamped neither side, intersecting
        p2 = np.array([1.0, 0.0, 0.0])
        q2 = np.array([0.0, 1.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 0.0, 10)

    def test_offset_vert(self):
        # clamped neither side, nonintersecting
        p2 = np.array([1.0, 0.0, 2.0])
        q2 = np.array([0.0, 1.0, 2.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 4.0, 10)

    def test_t(self):
        # these tests are all clamped on one side
        p2 = np.array([3.0, 0.0, 0.0])
        q2 = np.array([0.0, 3.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 0.5, 10)

        p2 = np.array([-2.0, 0.0, 0.0])
        q2 = np.array([0.0, -2.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 2.0, 10)

        p2 = np.array([1.0, 0.0, 0.0])
        q2 = np.array([2.0, -1.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 0.5, 10)

        p2 = np.array([0.0, 1.0, 0.0])
        q2 = np.array([-1.0, 2.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 0.5, 10)
    
    def test_parallel(self):
        # these tests are of exactly parallel lines, in reversed directions
        p2 = self.p1 + np.array([-1.0, 1.0, 0.0])
        q2 = self.q1 + np.array([-1.0, 1.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 2.0, 10)

        p2 = self.q1 + np.array([-1.0, 1.0, 0.0])
        q2 = self.p1 + np.array([-1.0, 1.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 2.0, 10)

    def test_degenerate_2(self):
        p2 = np.array([0.0, 1.0, 0.0])
        q2 = np.array([0.0, 1.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 0.5, 10)

    def test_degenerate_1(self):
        p2 = np.array([1.0, 0.0, 0.0])
        q2 = np.array([0.0, 1.0, 0.0])
        p1 = np.array([1.0, 1.0, 0.0])
        q1 = np.array([1.0, 1.0, 0.0])
        d2 = t.line_line(p1, q1, p2, q2)
        self.assertAlmostEqual(d2, 0.5, 10)

if __name__ == '__main__':
    unittest.main()
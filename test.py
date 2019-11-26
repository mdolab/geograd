import numpy as np 
import unittest
from triangles import triangles as t 

# normal case


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
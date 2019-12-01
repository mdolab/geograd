import numpy as np 
import unittest
from triangles import triangles as t
from triangles import triangles_db as tdb 
from triangles_complex import triangles as tcs 
h = 1e-15

def test_point_tri_cs(a, b, c, p, testcase):
    d2_base = t.point_tri(a, b, c, p)
    rev_derivs = tdb.point_tri_b(a, b, c, p, 1.0)

    # test derivs wrt a
    acs = a + np.array([(0+1j)*h, 0.0, 0.0])
    dd2_cs = np.imag(tcs.point_tri(acs, b, c, p)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[0][0], dd2_cs, 10)

    acs = a + np.array([0.0, (0+1j)*h, 0.0])
    dd2_cs = np.imag(tcs.point_tri(acs, b, c, p)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[0][1], dd2_cs, 10)

    acs = a + np.array([0.0, 0.0, (0+1j)*h])
    dd2_cs = np.imag(tcs.point_tri(acs, b, c, p)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[0][2], dd2_cs, 10)

    # test derivs wrt b
    bcs = b + np.array([(0+1j)*h, 0.0, 0.0])
    dd2_cs = np.imag(tcs.point_tri(a, bcs, c, p)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[1][0], dd2_cs, 10)

    bcs = b + np.array([0.0, (0+1j)*h, 0.0])
    dd2_cs = np.imag(tcs.point_tri(a, bcs, c, p)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[1][1], dd2_cs, 10)

    bcs = b + np.array([0.0, 0.0, (0+1j)*h])
    dd2_cs = np.imag(tcs.point_tri(a, bcs, c, p)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[1][2], dd2_cs, 10)

    # test derivs wrt c
    ccs = c + np.array([(0+1j)*h, 0.0, 0.0])
    dd2_cs = np.imag(tcs.point_tri(a, b, ccs, p)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[2][0], dd2_cs, 10)

    ccs = c + np.array([0.0, (0+1j)*h, 0.0])
    dd2_cs = np.imag(tcs.point_tri(a, b, ccs, p)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[2][1], dd2_cs, 10)

    ccs = c + np.array([0.0, 0.0, (0+1j)*h])
    dd2_cs = np.imag(tcs.point_tri(a, b, ccs, p)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[2][2], dd2_cs, 10)

    # test derivs wrt p
    pcs = p + np.array([(0+1j)*h, 0.0, 0.0])
    dd2_cs = np.imag(tcs.point_tri(a, b, c, pcs)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[3][0], dd2_cs, 10)

    pcs = p + np.array([0.0, (0+1j)*h, 0.0])
    dd2_cs = np.imag(tcs.point_tri(a, b, c, pcs)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[3][1], dd2_cs, 10)

    pcs = p + np.array([0.0, 0.0, (0+1j)*h])
    dd2_cs = np.imag(tcs.point_tri(a, b, c, pcs)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[3][2], dd2_cs, 10)

def test_line_line_cs(p1, q1, p2, q2, testcase):
    d2_base = t.line_line(p1, q1, p2, q2)
    rev_derivs = tdb.line_line_b(p1, q1, p2, q2, 1.0)

    # test derivs wrt p1
    p1cs = p1 + np.array([(0+1j)*h, 0.0, 0.0])
    dd2_cs = np.imag(tcs.line_line(p1cs, q1, p2, q2)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[0][0], dd2_cs, 10)

    p1cs = p1 + np.array([0.0, (0+1j)*h, 0.0])
    dd2_cs = np.imag(tcs.line_line(p1cs, q1, p2, q2)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[0][1], dd2_cs, 10)

    p1cs = p1 + np.array([0.0, 0.0, (0+1j)*h])
    dd2_cs = np.imag(tcs.line_line(p1cs, q1, p2, q2)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[0][2], dd2_cs, 10)

    # test derivs wrt q1
    q1cs = q1 + np.array([(0+1j)*h, 0.0, 0.0])
    dd2_cs = np.imag(tcs.line_line(p1, q1cs, p2, q2)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[1][0], dd2_cs, 10)

    q1cs = q1 + np.array([0.0, (0+1j)*h, 0.0])
    dd2_cs = np.imag(tcs.line_line(p1, q1cs, p2, q2)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[1][1], dd2_cs, 10)

    q1cs = q1 + np.array([0.0, 0.0, (0+1j)*h])
    dd2_cs = np.imag(tcs.line_line(p1, q1cs, p2, q2)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[1][2], dd2_cs, 10)

    # test derivs wrt p2
    p2cs = p2 + np.array([(0+1j)*h, 0.0, 0.0])
    dd2_cs = np.imag(tcs.line_line(p1, q1, p2cs, q2)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[2][0], dd2_cs, 10)

    p2cs = p2 + np.array([0.0, (0+1j)*h, 0.0])
    dd2_cs = np.imag(tcs.line_line(p1, q1, p2cs, q2)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[2][1], dd2_cs, 10)

    p2cs = p2 + np.array([0.0, 0.0, (0+1j)*h])
    dd2_cs = np.imag(tcs.line_line(p1, q1, p2cs, q2)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[2][2], dd2_cs, 10)

    # test derivs wrt q2
    q2cs = q2 + np.array([(0+1j)*h, 0.0, 0.0])
    dd2_cs = np.imag(tcs.line_line(p1, q1, p2, q2cs)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[3][0], dd2_cs, 10)

    q2cs = q2 + np.array([0.0, (0+1j)*h, 0.0])
    dd2_cs = np.imag(tcs.line_line(p1, q1, p2, q2cs)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[3][1], dd2_cs, 10)

    q2cs = q2 + np.array([0.0, 0.0, (0+1j)*h])
    dd2_cs = np.imag(tcs.line_line(p1, q1, p2, q2cs)-d2_base)/h
    # if (dd2_cs != 0.0): 
    #     print(dd2_cs)
    testcase.assertAlmostEqual(rev_derivs[3][2], dd2_cs, 10)

class TestPointTri(unittest.TestCase):
    def setUp(self):
        self.a = np.array([0.0, 0.0, 0.0])
        self.b = np.array([1.0, 0.0, 0.0])
        self.c = np.array([0.0, 1.0, 0.0])
    
    def test_interior_plane(self):
        p = np.array([0.25, 0.25, 0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.0, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

    def test_interior_nonplane(self):
        p = np.array([0.25, 0.25, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.25, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)


    def test_edges_plane(self):
        # edge ab
        p = np.array([0.5, -0.5, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.25, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # edge bc
        p = np.array([1.0, 1.0, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # edge ac
        p = np.array([-0.5, 0.5, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.25, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)


    def test_edges_nonplane(self):
        # edge ab
        p = np.array([0.5, -0.5, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # edge bc
        p = np.array([1.0, 1.0, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.75, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # edge ac
        p = np.array([-0.5, 0.5, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)


    def test_vertices_plane(self):
        # vertex a
        p = np.array([-0.5, -0.5, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # vertex b
        p = np.array([1.5, 0.0, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.25, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # vertex c
        p = np.array([0.0, 1.5, 0.0])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.25, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)


    def test_vertices_nonplane(self):
        # vertex a
        p = np.array([-0.5, -0.5, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.75, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # vertex b
        p = np.array([1.5, 0.0, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # vertex c
        p = np.array([0.0, 1.5, 0.5])
        d2 = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d2, 0.5, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

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
        test_line_line_cs(self.p1, self.q1, p2, q2, self)

    def test_offset_vert(self):
        # clamped neither side, nonintersecting
        p2 = np.array([1.0, 0.0, 2.0])
        q2 = np.array([0.0, 1.0, 2.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 4.0, 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


    def test_t(self):
        # these tests are all clamped on one side
        p2 = np.array([3.0, 0.0, 0.0])
        q2 = np.array([0.0, 3.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 0.5, 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


        p2 = np.array([-2.0, 0.0, 0.0])
        q2 = np.array([0.0, -2.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 2.0, 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


        p2 = np.array([1.0, 0.0, 0.0])
        q2 = np.array([2.0, -1.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 0.5, 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


        p2 = np.array([0.0, 1.0, 0.0])
        q2 = np.array([-1.0, 2.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 0.5, 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)

    
    def test_parallel(self):
        # these tests are of exactly parallel lines, in reversed directions
        p2 = self.p1 + np.array([-1.0, 1.0, 0.0])
        q2 = self.q1 + np.array([-1.0, 1.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 2.0, 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


        p2 = self.q1 + np.array([-1.0, 1.0, 0.0])
        q2 = self.p1 + np.array([-1.0, 1.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 2.0, 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


    def test_degenerate_2(self):
        p2 = np.array([0.0, 1.0, 0.0])
        q2 = np.array([0.0, 1.0, 0.0])
        d2 = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d2, 0.5, 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


    def test_degenerate_1(self):
        p2 = np.array([1.0, 0.0, 0.0])
        q2 = np.array([0.0, 1.0, 0.0])
        p1 = np.array([1.0, 1.0, 0.0])
        q1 = np.array([1.0, 1.0, 0.0])
        d2 = t.line_line(p1, q1, p2, q2)
        self.assertAlmostEqual(d2, 0.5, 10)
        test_line_line_cs(p1, q1, p2, q2, self)


if __name__ == '__main__':
    unittest.main()
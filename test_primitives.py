import numpy as np 
import unittest
from geograd import triangles as t
from geograd import triangles_db as tdb 
from geograd_complex import triangles as tcs 
h = 1e-15

def test_point_tri_cs(a, b, c, p, testcase):
    d_base = t.point_tri(a, b, c, p)
    rev_derivs = tdb.point_tri_b(a, b, c, p, 1.0)

    # test derivs wrt a
    acs = a + np.array([(0+1j)*h, 0.0, 0.0])
    dd_cs = np.imag(tcs.point_tri(acs, b, c, p)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[0][0], dd_cs, 10)

    acs = a + np.array([0.0, (0+1j)*h, 0.0])
    dd_cs = np.imag(tcs.point_tri(acs, b, c, p)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[0][1], dd_cs, 10)

    acs = a + np.array([0.0, 0.0, (0+1j)*h])
    dd_cs = np.imag(tcs.point_tri(acs, b, c, p)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[0][2], dd_cs, 10)

    # test derivs wrt b
    bcs = b + np.array([(0+1j)*h, 0.0, 0.0])
    dd_cs = np.imag(tcs.point_tri(a, bcs, c, p)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[1][0], dd_cs, 10)

    bcs = b + np.array([0.0, (0+1j)*h, 0.0])
    dd_cs = np.imag(tcs.point_tri(a, bcs, c, p)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[1][1], dd_cs, 10)

    bcs = b + np.array([0.0, 0.0, (0+1j)*h])
    dd_cs = np.imag(tcs.point_tri(a, bcs, c, p)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[1][2], dd_cs, 10)

    # test derivs wrt c
    ccs = c + np.array([(0+1j)*h, 0.0, 0.0])
    dd_cs = np.imag(tcs.point_tri(a, b, ccs, p)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[2][0], dd_cs, 10)

    ccs = c + np.array([0.0, (0+1j)*h, 0.0])
    dd_cs = np.imag(tcs.point_tri(a, b, ccs, p)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[2][1], dd_cs, 10)

    ccs = c + np.array([0.0, 0.0, (0+1j)*h])
    dd_cs = np.imag(tcs.point_tri(a, b, ccs, p)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[2][2], dd_cs, 10)

    # test derivs wrt p
    pcs = p + np.array([(0+1j)*h, 0.0, 0.0])
    dd_cs = np.imag(tcs.point_tri(a, b, c, pcs)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[3][0], dd_cs, 10)

    pcs = p + np.array([0.0, (0+1j)*h, 0.0])
    dd_cs = np.imag(tcs.point_tri(a, b, c, pcs)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[3][1], dd_cs, 10)

    pcs = p + np.array([0.0, 0.0, (0+1j)*h])
    dd_cs = np.imag(tcs.point_tri(a, b, c, pcs)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[3][2], dd_cs, 10)

def test_line_line_cs(p1, q1, p2, q2, testcase):
    d_base = t.line_line(p1, q1, p2, q2)
    rev_derivs = tdb.line_line_b(p1, q1, p2, q2, 1.0)

    # test derivs wrt p1
    p1cs = p1 + np.array([(0+1j)*h, 0.0, 0.0])
    dd_cs = np.imag(tcs.line_line(p1cs, q1, p2, q2)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[0][0], dd_cs, 10)

    p1cs = p1 + np.array([0.0, (0+1j)*h, 0.0])
    dd_cs = np.imag(tcs.line_line(p1cs, q1, p2, q2)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[0][1], dd_cs, 10)

    p1cs = p1 + np.array([0.0, 0.0, (0+1j)*h])
    dd_cs = np.imag(tcs.line_line(p1cs, q1, p2, q2)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[0][2], dd_cs, 10)

    # test derivs wrt q1
    q1cs = q1 + np.array([(0+1j)*h, 0.0, 0.0])
    dd_cs = np.imag(tcs.line_line(p1, q1cs, p2, q2)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[1][0], dd_cs, 10)

    q1cs = q1 + np.array([0.0, (0+1j)*h, 0.0])
    dd_cs = np.imag(tcs.line_line(p1, q1cs, p2, q2)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[1][1], dd_cs, 10)

    q1cs = q1 + np.array([0.0, 0.0, (0+1j)*h])
    dd_cs = np.imag(tcs.line_line(p1, q1cs, p2, q2)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[1][2], dd_cs, 10)

    # test derivs wrt p2
    p2cs = p2 + np.array([(0+1j)*h, 0.0, 0.0])
    dd_cs = np.imag(tcs.line_line(p1, q1, p2cs, q2)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[2][0], dd_cs, 10)

    p2cs = p2 + np.array([0.0, (0+1j)*h, 0.0])
    dd_cs = np.imag(tcs.line_line(p1, q1, p2cs, q2)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[2][1], dd_cs, 10)

    p2cs = p2 + np.array([0.0, 0.0, (0+1j)*h])
    dd_cs = np.imag(tcs.line_line(p1, q1, p2cs, q2)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[2][2], dd_cs, 10)

    # test derivs wrt q2
    q2cs = q2 + np.array([(0+1j)*h, 0.0, 0.0])
    dd_cs = np.imag(tcs.line_line(p1, q1, p2, q2cs)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[3][0], dd_cs, 10)

    q2cs = q2 + np.array([0.0, (0+1j)*h, 0.0])
    dd_cs = np.imag(tcs.line_line(p1, q1, p2, q2cs)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[3][1], dd_cs, 10)

    q2cs = q2 + np.array([0.0, 0.0, (0+1j)*h])
    dd_cs = np.imag(tcs.line_line(p1, q1, p2, q2cs)-d_base)/h
    # if (dd_cs != 0.0): 
    #     print(dd_cs)
    testcase.assertAlmostEqual(rev_derivs[3][2], dd_cs, 10)

def test_intersect_cs(a1, b1, c1, a2, b2, c2, testcase):
    d_base = t.intersect(a1, b1, c1, a2, b2, c2)
    rev_derivs = tdb.intersect_b(a1, b1, c1, a2, b2, c2, 1.0)

    # test derivs wrt a1
    a1cs = a1 + np.array([(0+1j)*h, 0.0, 0.0])
    d_cs = np.imag(tcs.intersect(a1cs, b1, c1, a2, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[0][0], d_cs, 10)

    a1cs = a1 + np.array([0.0, (0+1j)*h, 0.0])
    d_cs = np.imag(tcs.intersect(a1cs, b1, c1, a2, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[0][1], d_cs, 10)

    a1cs = a1 + np.array([0.0, 0.0, (0+1j)*h])
    d_cs = np.imag(tcs.intersect(a1cs, b1, c1, a2, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[0][2], d_cs, 10)

    # test derivs wrt b1
    b1cs = b1 + np.array([(0+1j)*h, 0.0, 0.0])
    d_cs = np.imag(tcs.intersect(a1, b1cs, c1, a2, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[1][0], d_cs, 10)

    b1cs = b1 + np.array([0.0, (0+1j)*h, 0.0])
    d_cs = np.imag(tcs.intersect(a1, b1cs, c1, a2, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[1][1], d_cs, 10)

    b1cs = b1 + np.array([0.0, 0.0, (0+1j)*h])
    d_cs = np.imag(tcs.intersect(a1, b1cs, c1, a2, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[1][2], d_cs, 10)

    # test derivs wrt c1
    c1cs = c1 + np.array([(0+1j)*h, 0.0, 0.0])
    d_cs = np.imag(tcs.intersect(a1, b1, c1cs, a2, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[2][0], d_cs, 10)

    c1cs = c1 + np.array([0.0, (0+1j)*h, 0.0])
    d_cs = np.imag(tcs.intersect(a1, b1, c1cs, a2, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[2][1], d_cs, 10)

    c1cs = c1 + np.array([0.0, 0.0, (0+1j)*h])
    d_cs = np.imag(tcs.intersect(a1, b1, c1cs, a2, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[2][2], d_cs, 10)

    # test derivs wrt a2
    a2cs = a2 + np.array([(0+1j)*h, 0.0, 0.0])
    d_cs = np.imag(tcs.intersect(a1, b1, c1, a2cs, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[3][0], d_cs, 10)

    a2cs = a2 + np.array([0.0, (0+1j)*h, 0.0])
    d_cs = np.imag(tcs.intersect(a1, b1, c1, a2cs, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[3][1], d_cs, 10)

    a2cs = a2 + np.array([0.0, 0.0, (0+1j)*h])
    d_cs = np.imag(tcs.intersect(a1, b1, c1, a2cs, b2, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[3][2], d_cs, 10)

    # test derivs wrt b2
    b2cs = b2 + np.array([(0+1j)*h, 0.0, 0.0])
    d_cs = np.imag(tcs.intersect(a1, b1, c1, a2, b2cs, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[4][0], d_cs, 10)

    b2cs = b2 + np.array([0.0, (0+1j)*h, 0.0])
    d_cs = np.imag(tcs.intersect(a1, b1, c1, a2, b2cs, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[4][1], d_cs, 10)

    b2cs = b2 + np.array([0.0, 0.0, (0+1j)*h])
    d_cs = np.imag(tcs.intersect(a1, b1, c1, a2, b2cs, c2)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[4][2], d_cs, 10)

    # test derivs wrt c2
    c2cs = c2 + np.array([(0+1j)*h, 0.0, 0.0])
    d_cs = np.imag(tcs.intersect(a1, b1, c1, a2, b2, c2cs)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[5][0], d_cs, 10)

    c2cs = c2 + np.array([0.0, (0+1j)*h, 0.0])
    d_cs = np.imag(tcs.intersect(a1, b1, c1, a2, b2, c2cs)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[5][1], d_cs, 10)

    c2cs = c2 + np.array([0.0, 0.0, (0+1j)*h])
    d_cs = np.imag(tcs.intersect(a1, b1, c1, a2, b2, c2cs)-d_base)/h
    # if (d_cs != 0.0): 
    #     print(d_cs)
    testcase.assertAlmostEqual(rev_derivs[5][2], d_cs, 10)

class TestPointTri(unittest.TestCase):
    def setUp(self):
        self.a = np.array([0.0, 0.0, 0.0])
        self.b = np.array([1.0, 0.0, 0.0])
        self.c = np.array([0.0, 1.0, 0.0])
    
    def test_interior_plane(self):
        p = np.array([0.25, 0.25, 0])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, 0.0, 10)
        #undefined derivative at exactly 0 dist due to the square root
        #test_point_tri_cs(self.a, self.b, self.c, p, self)

    def test_interior_nonplane(self):
        p = np.array([0.25, 0.25, 0.5])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, 0.5, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)


    def test_edges_plane(self):
        # edge ab
        p = np.array([0.5, -0.5, 0.0])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, 0.5, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # edge bc
        p = np.array([1.0, 1.0, 0.0])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, np.sqrt(0.5), 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # edge ac
        p = np.array([-0.5, 0.5, 0.0])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, 0.5, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)


    def test_edges_nonplane(self):
        # edge ab
        p = np.array([0.5, -0.5, 0.5])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, np.sqrt(0.5), 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # edge bc
        p = np.array([1.0, 1.0, 0.5])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, np.sqrt(0.75), 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # edge ac
        p = np.array([-0.5, 0.5, 0.5])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, np.sqrt(0.5), 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)


    def test_vertices_plane(self):
        # vertex a
        p = np.array([-0.5, -0.5, 0.0])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, np.sqrt(0.5), 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # vertex b
        p = np.array([1.5, 0.0, 0.0])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, 0.5, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # vertex c
        p = np.array([0.0, 1.5, 0.0])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, 0.5, 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)


    def test_vertices_nonplane(self):
        # vertex a
        p = np.array([-0.5, -0.5, 0.5])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, np.sqrt(0.75), 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # vertex b
        p = np.array([1.5, 0.0, 0.5])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, np.sqrt(0.5), 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

        # vertex c
        p = np.array([0.0, 1.5, 0.5])
        d = t.point_tri(self.a, self.b, self.c, p)
        self.assertAlmostEqual(d, np.sqrt(0.5), 10)
        test_point_tri_cs(self.a, self.b, self.c, p, self)

class TestLineLine(unittest.TestCase):
    def setUp(self):
        self.p1 = np.array([0.0,0.0,0.0])
        self.q1 = np.array([1.0, 1.0, 0.0])

    def test_intersect(self):
        # clamped neither side, intersecting
        p2 = np.array([1.0, 0.0, 0.0])
        q2 = np.array([0.0, 1.0, 0.0])
        d = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d, 0.0, 10)
        # cs derivatives are undefined when dist is exactly 0
        #test_line_line_cs(self.p1, self.q1, p2, q2, self)

    def test_offset_vert(self):
        # clamped neither side, nonintersecting
        p2 = np.array([1.0, 0.0, 2.0])
        q2 = np.array([0.0, 1.0, 2.0])
        d = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d, 2.0, 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


    def test_t(self):
        # these tests are all clamped on one side
        p2 = np.array([3.0, 0.0, 0.0])
        q2 = np.array([0.0, 3.0, 0.0])
        d = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d, np.sqrt(0.5), 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


        p2 = np.array([-2.0, 0.0, 0.0])
        q2 = np.array([0.0, -2.0, 0.0])
        d = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d, np.sqrt(2.0), 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


        p2 = np.array([1.0, 0.0, 0.0])
        q2 = np.array([2.0, -1.0, 0.0])
        d = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d, np.sqrt(0.5), 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


        p2 = np.array([0.0, 1.0, 0.0])
        q2 = np.array([-1.0, 2.0, 0.0])
        d = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d, np.sqrt(0.5), 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)

    
    def test_parallel(self):
        # these tests are of exactly parallel lines, in reversed directions
        p2 = self.p1 + np.array([-1.0, 1.0, 0.0])
        q2 = self.q1 + np.array([-1.0, 1.0, 0.0])
        d = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d, np.sqrt(2.0), 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


        p2 = self.q1 + np.array([-1.0, 1.0, 0.0])
        q2 = self.p1 + np.array([-1.0, 1.0, 0.0])
        d = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d, np.sqrt(2.0), 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


    def test_degenerate_2(self):
        p2 = np.array([0.0, 1.0, 0.0])
        q2 = np.array([0.0, 1.0, 0.0])
        d = t.line_line(self.p1, self.q1, p2, q2)
        self.assertAlmostEqual(d, np.sqrt(0.5), 10)
        test_line_line_cs(self.p1, self.q1, p2, q2, self)


    def test_degenerate_1(self):
        p2 = np.array([1.0, 0.0, 0.0])
        q2 = np.array([0.0, 1.0, 0.0])
        p1 = np.array([1.0, 1.0, 0.0])
        q1 = np.array([1.0, 1.0, 0.0])
        d = t.line_line(p1, q1, p2, q2)
        self.assertAlmostEqual(d, np.sqrt(0.5), 10)
        test_line_line_cs(p1, q1, p2, q2, self)

def test_intersect_permute_edges(a1, b1, c1, a2, b2, c2, answer, tol, testcase):
    d = t.intersect(a1, b1, c1, a2, b2, c2)
    testcase.assertAlmostEqual(d, answer, tol)
    test_intersect_cs(a1, b1, c1, a2, b2, c2, testcase)

    d = t.intersect(b1, c1, a1, a2, b2, c2)
    testcase.assertAlmostEqual(d, answer, tol)
    test_intersect_cs(b1, c1, a1, a2, b2, c2, testcase)

    d = t.intersect(c1, a1, b1, a2, b2, c2)
    testcase.assertAlmostEqual(d, answer, tol)
    test_intersect_cs(c1, a1, b1, a2, b2, c2, testcase)

    d = t.intersect(b1, a1, c1, a2, b2, c2)
    testcase.assertAlmostEqual(d, answer, tol)
    test_intersect_cs(b1, a1, c1, a2, b2, c2, testcase)

    d = t.intersect(a1, b1, c1, b2, c2, a2)
    testcase.assertAlmostEqual(d, answer, tol)
    test_intersect_cs(a1, b1, c1, b2, c2, a2, testcase)

    d = t.intersect(a1, b1, c1, c2, a2, b2)
    testcase.assertAlmostEqual(d, answer, tol)
    test_intersect_cs(a1, b1, c1, c2, a2, b2, testcase)

    d = t.intersect(a1, b1, c1, c2, b2, a2)
    testcase.assertAlmostEqual(d, answer, tol)
    test_intersect_cs(a1, b1, c1, c2, b2, a2, testcase)

class TestIntersect(unittest.TestCase):
    def setUp(self):
        self.a1 = np.array([0.0, 0.0, 0.0])
        self.b1 = np.array([1.0, 0.0, 0.0])
        self.c1 = np.array([0.0, 1.0, 0.0])


    def test_intersect_coplanar(self):
        a2 = self.a1 + np.array([2.0, 2.0, 0.0])
        b2 = np.array([1.0, 0.5, -0.1]) + np.array([2.0, 2.0, 0.0])
        c2 = np.array([0.5, 0.5, 5.0]) + np.array([2.0, 2.0, 0.0])
        test_intersect_permute_edges(self.a1, self.b1, self.c1, a2, b2, c2, 0.0, 10, self)

    def test_intersect_1_in_2(self):
        a2 = np.array([-1.0, 0.5, -0.1])
        b2 = np.array([1.0, 0.5, -0.1])
        c2 = np.array([0.5, 0.5, 5.0])
        test_intersect_permute_edges(self.a1, self.b1, self.c1, a2, b2, c2, 0.5, 10, self)

    def test_intersect_2_in_1(self):
        a2 = np.array([0.75, 0.5, -2.0])
        b2 = np.array([-0.25, 0.5, -2.0])
        c2 = np.array([0.25, 0.5, 1.0])
        test_intersect_permute_edges(self.a1, self.b1, self.c1, a2, b2, c2, 1/3, 10, self)

    def test_intersect_right_edge(self):
        a2 = np.array([1.0, 0.5, -2.0])
        b2 = np.array([0.0, 0.5, -2.0])
        c2 = np.array([0.5, 0.5, 1.0])
        test_intersect_permute_edges(self.a1, self.b1, self.c1, a2, b2, c2, 1/6, 10, self)

    def test_intersect_left_edge(self):
        a2 = np.array([0.5, 0.5, -2.0])
        b2 = np.array([-0.5, 0.5, -2.0])
        c2 = np.array([0.0, 0.5, 1.0])
        test_intersect_permute_edges(self.a1, self.b1, self.c1, a2, b2, c2, 1/6, 10, self)

    def test_intersect_scaling_independence(self):
        a2 = np.array([0.5, 0.5, -2.0])
        b2 = np.array([-0.5, 0.5, -2.0])
        c2 = np.array([0.0, 0.5, 1.0])
        test_intersect_permute_edges(self.a1*12, self.b1*10, self.c1*7, a2, b2, c2, 1/6, 10, self)

    def test_no_intersect_not_easy(self):
        a2 = np.array([3.0, 0.5, -0.1])
        b2 = np.array([5.0, 0.5, -0.1])
        c2 = np.array([4.5, 0.5, 5.0])
        test_intersect_permute_edges(self.a1, self.b1, self.c1, a2, b2, c2, 0.0, 10, self)

        a2 = np.array([-3.0, 0.5, -0.1])
        b2 = np.array([-5.0, 0.5, -0.1])
        c2 = np.array([-4.5, 0.5, 5.0])
        test_intersect_permute_edges(self.a1, self.b1, self.c1, a2, b2, c2, 0.0, 10, self)

    def test_no_intersect_easy(self):
        # separating axis test
        a2 = np.array([0.0, 0.0, 0.001])
        b2 = np.array([1.0, 0.0, 0.001])
        c2 = np.array([0.0, 1.0, 0.002])
        test_intersect_permute_edges(self.a1, self.b1, self.c1, a2, b2, c2, 0.0, 10, self)

        a2 = np.array([0.0, 0.0, -0.001])
        b2 = np.array([1.0, 0.0, -0.001])
        c2 = np.array([0.0, 1.0, -0.002])
        test_intersect_permute_edges(self.a1, self.b1, self.c1, a2, b2, c2, 0.0, 10, self)

if __name__ == '__main__':
    unittest.main()
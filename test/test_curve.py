import unittest
import sys
sys.path.insert(0, "isogenies-vdf")
import setup, curve, point, strategy
from sage.rings.integer_ring import Z as ZZ
from sage.misc.prandom import randint
from sage.arith.misc import next_prime
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class TestCurve(unittest.TestCase):

    def test_j_weierstrass(self):
        for s in setup.SETUPS.values():
            c = curve.Curve(s.alpha,s)
            self.assertEqual(c.j, c.weierstrass().j_invariant())
        
    def test_max_2_torsion(self):
        for s in setup.SETUPS.values():
            c = curve.Curve(s.alpha, s)
            n = c.max_2_torsion
            P = c.point_of_order(N=False, n=n, twist=False, deterministic=False)
            self.assertTrue(True) # ok if P is found
    
    def test_elligator(self, reps=10):
        for s in setup.SETUPS.values():
            c = curve.Curve(s.alpha, s)
            F = c.field
            for i in range(reps) :
                P = c.elligator(F.random_element(), twist=False)
                P_ws = P.weierstrass(c,twist=False)
                self.assertTrue(True) # ok if P_ws is created
    
    def test_elligator_twist(self, reps=10):
        for s in setup.SETUPS.values():
            c = curve.Curve(s.alpha, s)
            F = c.field
            self.assertTrue(F.is_prime_field())
            for i in range(reps) :
                P = c.elligator(F.random_element(), twist=True)
                P_ws = P.weierstrass(c,twist=True)
                self.assertTrue(True) # ok if P_ws is created
    
    #def test_point_of_order(self):
    #    # stuff with point_of_order(self, N=False, n=None, twist=False, deterministic=True)
    #    return True
        
    #def test__sqrt_expo(self):
    #    # stuff with _sqrt_expo(self, u)
    #    return True
    
    #def test__sqrt(self):
    #    # stuff with _sqrt(self, u, principal=True)
    #    return True
    
    def test_isogeny_forward(self, reps=10):
        s = setup.SETUPS['p14-toy']
        c = curve.Curve(s.alpha, s).to_gfp2()
        for i in range(reps):
            P = c.point_of_order(N=False, n=2, deterministic=False)
            while (2*P).x == 0 :
                P = c.point_of_order(N=False, n=2, deterministic=False)
            P2 = 2*P
            # P2 is a point of order 2 with non-zero x-coordinate
            alpha = P2.x/P2.z
            c.alpha = alpha
            c1, l = c.isogeny_forward((P,))
            self.assertEqual(l[0].weierstrass(c1).order(), 2)
            c1.alpha = l[0].x/l[0].z
            c2, l = c1.isogeny_forward(())
            self.assertEqual(c2.j, c.weierstrass().isogeny_codomain(P.weierstrass(c)).j_invariant())
            c = c2
    
    def test_isogeny_backward(self, reps=10):
        s = setup.SETUPS['p14-toy']
        c = curve.Curve(s.alpha, s).to_gfp2()
        P = c.point_of_order(N=False, n=1, deterministic=False)
        # P is a point of order 2
        while P.x == 0 :
            P = c.point_of_order(N=False, n=1, deterministic=False)
        self.assertTrue((2*P).is_zero() and not(P.is_zero()))
        for i in range(reps):
            Q = c.point_of_order(N=True,n=c.max_2_torsion,deterministic=False)
            c.alpha = P.x/P.z
            c1, l = c.isogeny_forward((Q,))
            ll = c.isogeny_backward(l[0])
            R = point.Point(ll[0].x,ll[0].z, c)
            self.assertEqual(R, 2*Q)
    
    def test_large_isogeny_forward(self):
        s = setup.SETUPS['p14-toy']
        c = curve.Curve(s.alpha, s).to_gfp2()
        for nn in range(1, c.max_2_torsion+1) :
            P = c.point_of_order(N=False, n=nn, deterministic=False)
            while ((2**(nn-1))*P).x == 0 :
                P = c.point_of_order(N=False, n=nn, deterministic=False)
            s = strategy.Strategy(0)
            c2, l = c.large_isogeny_forward(P, (), s.generate(nn-1, 1, 1))
            self.assertEqual(c2.j, c.weierstrass().isogeny_codomain(P.weierstrass(c)).j_invariant())

    #def test_large_isogeny_backward(self):
    #    s = setup.SETUPS['p14-toy']
    #    c = curve.Curve(s.alpha, s).to_gfp2()
    #    P = c.point_of_order(N=False, n=3, deterministic=False)
    #    c2, l = c.large_isogeny_forward(P, (), [1,1])        
    
if __name__ == '__main__':
    unittest.main()

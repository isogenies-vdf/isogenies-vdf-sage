import unittest
import sys
sys.path.insert(0, "isogenies-vdf")
import setup
import curve, point
from sage.rings.integer_ring import Z as ZZ
from sage.misc.prandom import randint
from sage.arith.misc import next_prime
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class TestCurve(unittest.TestCase):

    #def __init__(self, alpha, setup):
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
    
    def test_point_of_order(self):
        # stuff with point_of_order(self, N=False, n=None, twist=False, deterministic=True)
        # todo
        return True
        
    def test__sqrt_expo(self):
        # stuff with _sqrt_expo(self, u)
        return True
    
    def test__sqrt(self):
        # stuff wit _sqrt(self, u, principal=True)
        return True
        
    def test_isogeny_forward(self):
        # stuff with isogeny_forward(self, points, principal=True)
        return True
        
    def test_isogeny_backward(self):
        # stuff wit isogeny_backward(self, *points)
        return True
        
if __name__ == '__main__':
    unittest.main()

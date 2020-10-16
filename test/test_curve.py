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

    '''
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
        # stuff with _sqrt(self, u, principal=True)
        return True
    '''

    def test_isogeny_forward(self):
        return True
    
    def test_large_isogeny_forward(self):
        # stuff with isogeny_forward(self, points, principal=True)
        s = setup.SETUPS['p14-toy']
        c = curve.Curve(s.alpha, s).to_gfp2()
        P = c.point_of_order(N=False, n=2, deterministic=False)
        while ((2*P).x == 0) :
            P = c.point_of_order(N=False, n=2, deterministic=False)
        # P is a point of order 4
        assert (4*P).is_zero()
        assert not((2*P).is_zero())
        assert not(P.is_zero())
        P2 = 2*P
        # P2 is a point of order 2
        assert (2*P2).is_zero()
        assert not(P2.is_zero())
        j1 = c.weierstrass().isogeny_codomain(P.weierstrass(c)).j_invariant()

        c2, l = c.isogeny_forward((), principal=False)
        c2, l = c2.isogeny_forward((), principal=False)
        while(c2.weierstrass().j_invariant() != j1) :
            c2, l = c.isogeny_forward((), principal=True)
            c2, l = c2.isogeny_forward((), principal=True)
            print(c2.weierstrass().j_invariant(), j1)
        print('ok')
        #NOT WORKING
        '''
        E1, l1 = c.large_isogeny_forward(P, (), [1], 2, principal=False)
        P2 = 2*P
        print('ker=', P2.x/P2.z)
        E2, l2 = c.isogeny_forward((P,), principal=False, alpha = P2.x/P2.z)
        fP = l2[0]
        print(E2.weierstrass().j_invariant())
        print('fP=', fP.weierstrass(E2))
        print('order=', fP.weierstrass(E2).order().factor())
        #E4, l4 = E2.isogeny_forward((), principal=False, alpha=fP.x/fP.z)
        #this is not working
        #print(E1.j)
        #print(E4.j)
        '''
        return True
        
    def test_isogeny_backward(self):
        # stuff wit isogeny_backward(self, *points)
        return True
        
if __name__ == '__main__':
    unittest.main()

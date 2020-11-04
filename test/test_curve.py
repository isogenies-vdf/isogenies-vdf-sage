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

    #def test_j_weierstrass(self):
    #    for s in setup.SETUPS.values():
    #        c = curve.Curve(s.alpha,s)
    #        self.assertEqual(c.j, c.weierstrass().j_invariant())
    #    
    #def test_max_2_torsion(self):
    #    for s in setup.SETUPS.values():
    #        c = curve.Curve(s.alpha, s)
    #        n = c.max_2_torsion
    #        P = c.point_of_order(N=False, n=n, twist=False, deterministic=False)
    #        self.assertTrue(True) # ok if P is found
    #
    #def test_elligator(self, reps=10):
    #    for s in setup.SETUPS.values():
    #        c = curve.Curve(s.alpha, s)
    #        F = c.field
    #        for i in range(reps) :
    #            P = c.elligator(F.random_element(), twist=False)
    #            P_ws = P.weierstrass(c,twist=False)
    #            self.assertTrue(True) # ok if P_ws is created
    #
    #def test_elligator_twist(self, reps=10):
    #    for s in setup.SETUPS.values():
    #        c = curve.Curve(s.alpha, s)
    #        F = c.field
    #        self.assertTrue(F.is_prime_field())
    #        for i in range(reps) :
    #            P = c.elligator(F.random_element(), twist=True)
    #            xx = P.x/P.z
    #            P_ws = P.weierstrass(c,twist=True)
    #            self.assertTrue(True) # ok if P_ws is created

    
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
        c = curve.Curve(s.alpha, s)
        for c in [c, c.to_gfp2()]:
            for i in range(reps):
                P = c.point_of_order(N=False, n=1, deterministic=False)
                while P.x == 0 :
                    P = c.point_of_order(N=False, n=1, deterministic=False)
                alpha = P.x/P.z
                c.alpha = alpha
                R = c.point_of_order(N=True)
                c1, l = c.isogeny_forward((R,))
                Pws = P.weierstrass(c)
                phi = c.weierstrass().isogeny(Pws)
                self.assertEqual(phi.codomain().j_invariant(), c1.j)
                c = c1

    def test_isogeny_backward(self, reps=10):
        '''
        Checks that the composition of forward and backward is the multiplication by 2
        '''
        s = setup.SETUPS['p14-toy']
        c = curve.Curve(s.alpha, s)
        for c in [c, c.to_gfp2()]:
            for j in range(reps):
                P = c.point_of_order(N=False, n=1, deterministic=False)
                # P is a point of order 2
                while P.x == 0 :
                    P = c.point_of_order(N=False, n=1, deterministic=False)
                self.assertTrue((2*P).is_zero() and not(P.is_zero()))
                c.alpha = P.x/P.z
                Q = c.point_of_order(N=True,n=c.max_2_torsion,deterministic=False)
                c1, l = c.isogeny_forward((Q,))
                ll = c.isogeny_backward((l[0],))
                R = point.Point(ll[0].x,ll[0].z, c)
                self.assertEqual(R, 2*Q)
                c = c1
    
    def test_large_isogeny_forward(self):
        s = setup.SETUPS['p14-toy']
        c = curve.Curve(s.alpha, s)
        for c in [c, c.to_gfp2()]:
            for nn in range(2, c.max_2_torsion+1) :
                R = c.point_of_order(N=True,n=0,deterministic=False)
                s = strategy.Strategy(0)
                c2, l, P = c.large_isogeny_forward((R,), s.generate(nn, 1, 1))
                self.assertEqual(c2.j, c.weierstrass().isogeny_codomain(P.weierstrass(c)).j_invariant())
                                
    def test_large_isogeny_backward_with_dual(self, reps = 10):
        s = setup.SETUPS['p14-toy']
        c = curve.Curve(s.alpha, s)
        for c in [c, c.to_gfp2()]:
            for nn in range(2, c.max_2_torsion+1) :
                R = c.point_of_order(N=True,n=0,deterministic=False)
                s = strategy.Strategy(0)
                c2, l, P, list_of_curves = c.large_isogeny_forward((R,), s.generate(nn, 1, 1), with_dual=True)
                KernelDual = l[0]
                self.assertTrue(((2**nn) * KernelDual).is_zero())
                self.assertFalse(((2**(nn-1)) * KernelDual).is_zero())
                phiR = l[1]
                ll = c2.large_isogeny_backward(KernelDual, (phiR,), s.generate(nn,1,1), list_of_curves)
                hatphi_phi_R = ll[0]
                self.assertEqual(hatphi_phi_R, (2**nn)*R)

if __name__ == '__main__':
    unittest.main()

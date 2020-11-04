import unittest
import sys
sys.path.insert(0, "isogenies-vdf")
from setup import TrustedSetup
import tate, point
from sage.rings.integer_ring import Z as ZZ
from sage.misc.prandom import randint

class TestTate(unittest.TestCase):

    def test_double_line_jac(self, reps=10):
        s = TrustedSetup(f=1, n=8, N=53, alpha=5052)
        E = s.E0.to_gfp2().weierstrass()
        a2 = E.a2()
        pp = tate.Tate(0)
        for i in range(reps) :
            S = E.random_point()
            P = E.random_point()
            l, doublePoint = pp.double_line_jac(S, P, a2)
            if doublePoint[2] == 0 :
                self.assertEqual((2*S)[2], 0)
            else :
                #normalization of the point
                doublePoint = (doublePoint[0]/doublePoint[2]**2, doublePoint[1]/doublePoint[2]**3)
                self.assertEqual(doublePoint[0], (2*S)[0])
                self.assertEqual(doublePoint[1], (2*S)[1])
            self.assertEqual(l[0]/l[1], S._line_(S, P))
        
    def test_add_line_jac(self, reps=10):
        s = TrustedSetup(f=1, n=8, N=53, alpha=5052)
        E = s.E0.to_gfp2().weierstrass()
        a2 = E.a2()
        pp = tate.Tate(0)
        for i in range(reps) :
            S = E.random_point()
            P = E.random_point()
            Q = E.random_point()
            l, addPoint = pp.add_line_jac(S, P, Q, a2)
            if addPoint[2] == 0 :
                self.assertEqual((S+P)[2], 0)
            else :
                #normalization of the point
                addPoint = (addPoint[0]/addPoint[2]**2, addPoint[1]/addPoint[2]**3)
                self.assertEqual(addPoint[0], (S+P)[0])
                self.assertEqual(addPoint[1], (S+P)[1])
            self.assertEqual(l[0]/l[1], S._line_(P, Q))
            
    def test_vertical_line_jac(self, reps=10):
        s = TrustedSetup(f=1, n=8, N=53, alpha=5052)
        E = s.E0.to_gfp2().weierstrass()
        a2 = E.a2()
        pp = tate.Tate(0)
        for i in range(reps) :
            S = E.random_point()
            Q = E.random_point()
            v = pp.vertical_line_jac(S, Q)
            self.assertEqual(v[0]/v[1], S._line_(-S, Q))

    def test_miller(self, reps=10) :
        s = TrustedSetup(f=1, n=8, N=53, alpha=5052)
        E = s.E0.to_gfp2().weierstrass()
        a2 = E.a2()
        pp = tate.Tate(0)
        for i in range(100) :
            S = E.random_point()
            Q = E.random_point()
            millerLoop = ZZ(randint(1,s.N-1))
            mil, P = pp.miller(S, Q, millerLoop, a2)
            if P[2] == 0 :
                self.assertEqual((millerLoop*S)[2], 0)
            else :
                #normalization of the point
                P = (P[0]/P[2]**2, P[1]/P[2]**3)
                self.assertEqual(P[0], (millerLoop*S)[0])
                self.assertEqual(P[1], (millerLoop*S)[1])
            self.assertEqual(s.GFp2(mil[0]/mil[1]), S._miller_(Q, millerLoop))
        
    def test_tate_denominator(self, reps=10) :
        s = TrustedSetup(f=1, n=8, N=53, alpha=5052)
        E = s.E0.to_gfp2().weierstrass()
        a2 = E.a2()
        pp = tate.Tate(0)
        for i in range(reps) :
            P = ((s.p+1)//s.N) * E.random_point()
            while(P.order() != s.N) :
                P = ((s.p+1)//s.N) * E.random_point()
            Q = ((s.p+1)//s.N) * E.random_point()
            while(Q.order() != s.N or P.weil_pairing(Q, s.N) == 1) :
                Q = ((s.p+1)//s.N) * E.random_point()
            t = pp.tate(P, Q, s)
            self.assertEqual(t, P.tate_pairing(Q, s.N, 2))

    def test_tate_no_denominator(self, reps=10):
        s = TrustedSetup(f=1, n=8, N=53, alpha=5052)
        E_p = s.E0.weierstrass()
        E_p2 = s.E0.to_gfp2().weierstrass()
        a2 = E_p2.a2()
        pp = tate.Tate(0)
        for i in range(reps) :
            P = ((s.p+1)//s.N) * E_p.random_point()
            while(P.order() != s.N) :
                P = ((s.p+1)//s.N) * E_p.random_point()
            P = E_p2(P)
            Q = s.E0.point_of_order(N=True, n=0, twist=True, deterministic=False).weierstrass(s.E0, twist=True)
            ii = E_p2.base_field().gen()
            QQ = E_p2([-Q[0], ii*Q[1]])
            Q = QQ
            t = pp.tate(P, Q, s, denominator=False)
            self.assertEqual(t, P.tate_pairing(Q, s.N, 2))
            
if __name__ == '__main__':
    unittest.main()

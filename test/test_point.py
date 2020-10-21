import unittest
import sys
sys.path.insert(0, "isogenies-vdf")
import curve, point, setup
from sage.rings.integer_ring import Z as ZZ
from sage.misc.prandom import randint

class TestPoint(unittest.TestCase):

    def test_normalize(self, reps=5) :
        for s in setup.SETUPS.values():
            c = curve.Curve(s.alpha,s)
            IE = c.iterelligator(deterministic=False)
            for i in range(reps):
                P = next(IE)
                Pn = P.normalize()
                self.assertEqual(P, Pn)
                    
    def test_dbl(self, reps=5):
        for s in setup.SETUPS.values():
            c = curve.Curve(s.alpha,s)
            IE = c.iterelligator(deterministic=False)
            for i in range(reps):
                P = next(IE)
                P_ws = P.weierstrass(c)
                DoubleP = 2*P
                self.assertTrue((2*P_ws)[0]/(2*P_ws)[2] == DoubleP.x/DoubleP.z)
                    
    #def test_add(self) :
    #    # ?
    #    return True

    def test_mul(self, reps=5):
        for s in setup.SETUPS.values():
            c = curve.Curve(s.alpha,s)
            IE = c.iterelligator(deterministic=False)
            for i in range(reps):
                P = next(IE)
                P_ws = P.weierstrass(c)
                k = randint(0, 2**s.n * s.N * s.f - 2)
                kP = k*P
                self.assertTrue((k*P_ws)[0]/(k*P_ws)[2] == kP.x/kP.z)

    def test_trace(self, reps=5):
        for s in setup.SETUPS.values():
            c = curve.Curve(s.alpha,s)
            IE = c.iterelligator(deterministic=False)
            for i in range(reps):
                P = next(IE)
                tP = P.trace(c)
                P_ws = P.weierstrass(c)
                piP_ws = P_ws.curve()(P_ws[0]**s.p, P_ws[1]**s.p, P_ws[2]**s.p)
                tP_ws = P_ws + piP_ws
                self.assertEqual(tP.x/tP.z, tP_ws[0]/tP_ws[2])
                    
    #def test_weierstrass(self):
    #    #?
    #    return True
    
if __name__ == '__main__':
    unittest.main()

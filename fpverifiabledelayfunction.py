# -*- coding: utf-8 -*- 
from point import Point
from sage.rings.integer_ring import ZZ
from copy import copy
import pairing
from verifiabledelayfunction import VerifiableDelayFunction

class FpVerifiableDelayFunction(VerifiableDelayFunction):
    def __init__(self, method, strategy, curve, delay):
        #change it with super().__init__(method, strategy, curve, delay)
        # For the moment, python3 does not work...
        self.method = method
        self.strategy = strategy
        self.curve = curve
        Delta = delay
        while Delta % (curve.n-2) != 0 or (Delta // (curve.n-2)) % 2 != 0:
            Delta += 1
        self.delay = Delta

    def setup(self) :
        '''
        OUTPUT:
        * P the first point of the protocol
        * the ouput of isogeny_walk
        '''

        # point P of order N
        P = self.curve.cof_P * self.curve.random_point(1, True)
        while P.z == 0 :
            P = self.curve.cof_P * self.curve.random_point(1, True)
        return [P] + P.isogeny_walk(self.delay, self.strategy, self.method, 'fp')

    def evaluate(self, Q, curvesPath, kernelsOfBigSteps):
        '''
        INPUT:
        * Q the second point of the protocol
        * curvesPath from the setup
        * kernelsOfBigSteps from the setup
        OUTPUT:
        * hat_phiQ the image of Q by the dual walk
        '''

        if self.delay % (self.curve.n-2) != 0 :
            raise RuntimeError('Delta is not a multiple of n-2')
        k = ZZ((self.curve.n-2)//2)
        T = Q
        c_t = copy(Q.curve)

        if self.method == 'kernel4k' :
            for R in kernelsOfBigSteps :
                R = R.change_iso_curve(c_t.a)
                [T, kernelPoint, listOfCurves] = R.isogeny_degree4k(T, k, 'withoutKernel', self.strategy, stop=1)
                c_t = T.curve
        elif self.method == 'kernel4' :
            cpt = 0
            for c1 in curvesPath :
                R = Point(1,1,c1).change_iso_curve(c_t.a)
                [T] = R.isogeny_degree4([T])
                c_t = T.curve
            cpt += 1
        return T.change_iso_curve(self.curve.a)

    def verify(self, P, phiP, Q, hat_phiQ) :
        '''
        INPUT:
        * P the first point of the protocol
        * phiP the image of P
        * Q the second point of the protocol
        * hat_phiQ dual image of Q
        OUTPUT:
        * true/false depending on the verification
        '''

        if not(hat_phiQ.in_curve() and hat_phiQ.is_prime_order_point(self.curve.N)) :
            raise RuntimeError('evaluation step does not give point of the curve of order N')

        # this does not depend on the eval answer, and can be computed before the eval
        P_ws = P.weierstrass()
        phiP_ws = phiP.weierstrass()
        #this needs to be computed here
        Q_ws = Q.weierstrass()
        hat_phiQ_ws = hat_phiQ.weierstrass()

        _Z, mil11 = pairing.miller(hat_phiQ_ws, P_ws, ZZ(self.curve.N), denominator=False)
        e1 = pairing.exponentiation(self.curve, mil11[0]/mil11[1])

        _Z, mil22 = pairing.miller(Q_ws, phiP_ws, ZZ(self.curve.N), denominator=False)
        e2 = pairing.exponentiation(self.curve, mil22[0]/mil22[1])

        if e1 != 1 :
            if e1 == e2 :
                return True
            if e1 == 1/e2:
                return True
            # Pairing equation does not hold
            return False
        # e_Tr_hat_phiQ_P = 1
        return False

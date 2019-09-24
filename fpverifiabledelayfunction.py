# -*- coding: utf-8 -*- 
from sage.rings.integer_ring import ZZ
from copy import copy
import curve
from point import Point
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
        P = self.curve.pairing_group_random_point(extension_degree=1, twist=True)
        return [P] + self.setup_walk(1, P, stop=1)

    def evaluate(self, Q, dualKernels):
        '''
        INPUT:
        * Q the second point of the protocol
        * dualKernels from the setup
        OUTPUT:
        * hat_phiQ the image of Q by the dual walk
        '''
        return self.evaluation_walk(Q, dualKernels, 1)

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

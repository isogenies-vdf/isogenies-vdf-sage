# -*- coding: utf-8 -*- 
import logging
from sage.rings.integer_ring import ZZ
from copy import copy
import curve
from point import Point
import pairing
from verifiabledelayfunction import VerifiableDelayFunction

class FpVerifiableDelayFunction(VerifiableDelayFunction):
    def __init__(self, method, curve, delay):
        #change it with super().__init__(method, strategy, curve, delay)
        # For the moment, python3 does not work...
        self.method = method
        self.curve = curve
        Delta = delay
        while Delta % (curve.n-2) != 0 or (Delta // (curve.n-2)) % 2 != 0:
            Delta += 1
        self.delay = Delta
        # this is the strategy computation given in Luca De Feo's answer on crypto.stackexchange.com
        # it could be hard-coded.
        S = { 1: [] }
        C = { 1: 0 }
        p = 1
        q = 1
        nbIsog = (curve.n-2)//2
        for i in range(2, nbIsog+2):
            b, cost = min(((b, C[i-b] + C[b] + b*p + (i-b)*q) for b in range(1,i)), key=lambda t: t[1])
            S[i] = [b] + S[i-b] + S[b]
            C[i] = cost
        self.strategy = S[nbIsog+1]

    def setup(self) :
        P = self.curve.pairing_group_random_point(extension_degree=1, twist=True)
        logging.debug('P = %s', str(P))
        return [P] + self.setup_walk(1, P, stop=1, conditionP4=False)

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

        logging.debug('No denominator computed.')

        _Z, mil11 = pairing.miller(hat_phiQ_ws, P_ws, ZZ(self.curve.N), denominator=False)
        e1 = pairing.exponentiation(self.curve, mil11[0]/mil11[1])
        logging.debug('f_{N, hatphiQ}(P) = %s', str(mil11))
        logging.debug('e(hatphiQ, P) = %s', str(e1))


        _Z, mil22 = pairing.miller(Q_ws, phiP_ws, ZZ(self.curve.N), denominator=False)
        e2 = pairing.exponentiation(self.curve, mil22[0]/mil22[1])
        logging.debug('f_{N, Q}(phiP) = %s', str(mil22))
        logging.debug('e(Q, phiP) = %s', str(e2))

        if e1 != 1 :
            if e1 == e2 :
                return True
            if e1 == 1/e2:
                return True
            # Pairing equation does not hold
            logging.debug('Pairing equation does not hold.')
            return False
        # e_hat_phiQ_P = 1
        logging.debug('e(HatPhiQ, P) = 1.')
        return False

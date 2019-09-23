# -*- coding: utf-8 -*- 
import pairing
from sage.rings.integer_ring import ZZ
from point import Point
import curve
from copy import copy
from verifiabledelayfunction import VerifiableDelayFunction

class Fp2VerifiableDelayFunction(VerifiableDelayFunction):
    def __init__(self, method, strategy, curve, delay):
        #change it with super().__init__(method, strategy, curve, delay)
        # For the moment, python3 does not work...
        self.method = method
        self.strategy = strategy
        self.curve = curve
        Delta = delay
        while Delta % curve.n != 0 or (Delta // curve.n) % 2 != 0:
            Delta += 1
        self.delay = Delta

    def setup(self):
        P = self.curve.pairing_group_random_point(extension_degree=1, twist=True)
        return [P] + P.random_isogeny_walk(self.delay, 2, self.strategy, 0, self.method)

    def evaluate(self, Q, curvesPath, kernelsOfBigSteps):
        '''
        INPUT:
        * Q the second point of the protocol
        * curvesPath from the setup
        * kernelsOfBigSteps from the setup
        OUTPUT:
        * Tr_hat_phiQ the list of the possible images of Q by the dual walk composed by Trace (4 possible because of sign pb)
        '''
        if self.delay % self.curve.n != 0 :
            raise RuntimeError('Delta is not a multiple of n')
        T = Q.isogeny_walk(curvesPath, kernelsOfBigSteps, self.strategy)
        #T = hatphi(Q)

        #Trace trick
        frob_T = Point(T.x**c.p, T.z**c.p, T.curve)

        #not efficient
        fQ_ws = frob_T.weierstrass()
        Q_ws = T.weierstrass()

        R1 = Q_ws + fQ_ws
        R2 = Q_ws - fQ_ws

        #the (+/-) point to return is the one defined over Fp :-)
        return self.curve.getPointFromWeierstrass(R1) if R1[0] in self.curve.Fp and R1[1] in self.curve.Fp else self.curve.getPointFromWeierstrass(R2)

    def verify(self, P, phiP, Q, Tr_hat_phiQ) :
        '''
        INPUT:
        * P the first point of the protocol
        * phiP the image of P
        * Q the second point of the protocol
        * Tr_hat_phiQ the list of hat_phiQ + frob(hat_phiQ) and hat_phiQ - frob(hat_phiQ)
        OUTPUT:
        * true/false depending on the verification
        '''

        if not(Tr_hat_phiQ.in_curve() and Tr_hat_phiQ.x in c.Fp and Tr_hat_phiQ.z in c.Fp) :
            raise RuntimeError('evaluation step does not give point of the curve defined over Fp')

        # this does not depend on the eval answer, can be computed before the eval
        P_ws = P.weierstrass()
        phiP_ws = phiP.weierstrass()
        #this needs to be computed here
        Q_ws = Q.weierstrass()
        Tr_hat_phiQ_ws = Tr_hat_phiQ.weierstrass()

        _Z, mil11 = pairing.miller(Tr_hat_phiQ_ws, P_ws, ZZ(c.N), denominator=True)
        e1 = pairing.exponentiation(c, mil11[0]/mil11[1])

        _Z, mil22 = pairing.miller(Q_ws, phiP_ws, ZZ(c.N), denominator=True)
        e2_squared = pairing.exponentiation(c, mil22[0]/mil22[1])**2

        if e1 != 1 :
            if e1 == e2_squared :
                return True
            if e1 == 1/e2_squared:
                return True
            print('pairing eq doesnt hold')
            # Pairing equation does not hold
            return False
        print('pairing eq 1')
        # e_Tr_hat_phiQ_P = 1
        return False

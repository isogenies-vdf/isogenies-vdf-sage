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

    def setup(self) : #ADD ARGUMENTS!
        '''
        INPUT:
        * curve the elliptic curve
        * verbose for the comments
        * method for the method of storing the isogeny (kernel4, kernel4k or withoutKernel)
        OUTPUT:
        * P the first point of the protocol
        * the ouput of isogeny_walk
        '''

        # point P of order N
        P = curve.cof_P * curve.random_point(1, True)
        while P.z == 0 :
            P = curve.cof_P * curve.random_point(1, True)
        return [P] + isogeny_walk(curve, P, verbose, method)

    def evaluate(self, c, setup, Q, verbose, method): #CHANGE ARGUMENTS!
        '''
        INPUT:
        * c the elliptic curve
        * setup the setup from vdf_setup function
        * Q the second point of the protocol
        * verbose for the comments
        * method for the method of storing the isogeny (kernel4, kernel4k or withoutKernel)
        OUTPUT:
        * Tr_hat_phiQ the list of the possible images of Q by the dual walk composed by Trace (4 possible because of sign pb)
        '''
        [P, c_prime, curvesPath, kernelsOfBigSteps, phiP] = setup
        if c.Delta % c.n != 0 :
            raise RuntimeError('Delta is not a multiple of n')
        k = ZZ(c.n//2)

        T = Q
        # At the moment, isogeny codomain is "up to isomorphism".
        # From the kernels given in setup, I need to move them into the right curve.
        # That is why `kernel4` is not efficient (change_iso_curve times the number
        # of steps in the walk on the graph. In kernel4k, it is reduced by a factor
        # 1000 .
        if method == 'kernel4k' :
            for R in kernelsOfBigSteps :
                R = R.change_iso_curve(T.curve.a)
                [T, kernelPoint, listOfCurves] = R.isogeny_degree4k(T, k, method='withoutKernel', strategy=self.strategy)
        elif method == 'kernel4' :
            for c1 in curvesPath:
                R = Point(1, 1, c1).change_iso_curve(T.curve.a)
                [T] = R.isogeny_degree4([T])

        T = T.change_iso_curve(c.a)
        #T = hatphi(Q)

        #Trace trick
        frob_T = Point(T.x**c.p, T.z**c.p, T.curve)

        #not efficient
        fQ_ws = frob_T.weierstrass()
        Q_ws = T.weierstrass()

        R1 = Q_ws + fQ_ws
        R2 = Q_ws - fQ_ws

        #the (+/-) point to return is the one defined over Fp :-)
        return c.getPointFromWeierstrass(R1) if R1[0] in c.Fp and R1[1] in c.Fp else c.getPointFromWeierstrass(R2)

    def verify(self, c, setup, Q, Tr_hat_phiQ) :
        '''
        INPUT:
        * c the elliptic curve
        * setup the setup from the vdf_setup function
        * Q the second point of the protocol
        * Tr_hat_phiQ the list of hat_phiQ + frob(hat_phiQ) and hat_phiQ - frob(hat_phiQ)
        OUTPUT:
        * true/false depending on the verification
        '''
        [P, c_prime, curvesPath, kernelsOfBigSteps, phiP] = setup

        if not(Tr_hat_phiQ.in_curve() and Tr_hat_phiQ.x in c.Fp and Tr_hat_phiQ.z in c.Fp) :
            raise RuntimeError('evaluation step does not give point of the curve defined over Fp')

        # this does not depend on the eval answer, can be computed before the eval
        P_ws = P.weierstrass()
        phiP_ws = phiP.weierstrass()
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
            # Pairing equation does not hold
            return False
        # e_Tr_hat_phiQ_P = 1
        return False

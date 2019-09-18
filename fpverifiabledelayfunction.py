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

    def setup(self, curve, verbose, method) :
        '''
        INPUT:
        * curve the initial elliptic curve
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

    def evaluate(self, c, setup, Q, verbose, method):
        '''
        INPUT:
        * c the initial elliptic curve
        * setup the setup from vdf_setup function
        * Q the second point of the protocol
        * verbose
        * method for the method of storing the isogeny (kernel4, kernel4k or withoutKernel)
        OUTPUT:
        * hat_phiQ the image of Q by the dual walk
        '''
        [P, c_prime, curvesPath, kernelsOfBigSteps, phiP] = setup
        if c.Delta % (c.n-2) != 0 :
            raise RuntimeError('Delta is not a multiple of n-2')
        k = ZZ((c.n-2)//2)
        T = Q
        c_t = copy(c_prime)

        if method == 'kernel4k' :
            for R in kernelsOfBigSteps :
                R = R.change_iso_curve(c_t.a)
                [T, kernelPoint, listOfCurves] = R.isogeny_degree4k(T, k, method='withoutKernel', strategy=self.strategy, stop=1)
                c_t = T.curve
        elif method == 'kernel4' :
            cpt = 0
            for c1 in curvesPath :
                R = Point(1,1,c1).change_iso_curve(c_t.a)
                [T] = R.isogeny_degree4([T])
                c_t = T.curve
            cpt += 1
        return T.change_iso_curve(c.a)

    def verify(self, curve, setup, Q, hat_phiQ) :
        '''
        INPUT:
        * curve the elliptic curve
        * setup the setup from the vdf_setup function
        * Q the second point of the protocol
        * hat_phiQ dual image of Q
        OUTPUT:
        * true/false depending on the verification
        '''

        #we just need P and phiP for the moment... [P, curve_prime, curvesPath, kernelsOfBigSteps, phiP] = setup
        P, phiP = setup [0], setup[-1]

        if not(hat_phiQ.in_curve() and hat_phiQ.is_prime_order_point(curve.N)) :
            raise RuntimeError('evaluation step does not give point of the curve of order N')

        # this does not depend on the eval answer, and can be computed before the eval
        P_ws = P.weierstrass()
        phiP_ws = phiP.weierstrass()
        #this needs to be computed here
        Q_ws = Q.weierstrass()
        hat_phiQ_ws = hat_phiQ.weierstrass()

        _Z, mil11 = pairing.miller(hat_phiQ_ws, P_ws, ZZ(curve.N), denominator=False)
        e1 = pairing.exponentiation(curve, mil11[0]/mil11[1])

        _Z, mil22 = pairing.miller(Q_ws, phiP_ws, ZZ(curve.N), denominator=False)
        e2 = pairing.exponentiation(curve, mil22[0]/mil22[1])

        if e1 != 1 :
            if e1 == e2 :
                return True
            if e1 == 1/e2:
                return True
            # Pairing equation does not hold
            return False
        # e_Tr_hat_phiQ_P = 1
        return False

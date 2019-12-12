# -*- coding: utf-8 -*-
import curve, point, pairing
import logging

class VerifiableDelayFunction():
    '''
    A superclass for all variants of the VDF.
    '''
    
    def __init__(self, setup, delay):
        '''
        Initialize VDF from given setup, with given delay.
        '''
        self.setup = setup
        self.delay = delay
        
        self.P = self.setup.E0.point_of_order(N=True, n=0, twist=True, deterministic=True)
        self.E1, self.fP, self.walk_data = self._setup_walk()

    def random_input(self):
        'Sample a random input point for the VDF'
        return self.E1.point_of_order(N=True, n=0, twist=False, deterministic=False)
        
    def __repr__(self):
        return 'Verifiable delay function with delay %d' % self.delay

    __str__ = __repr__


class VDF_GFp(VerifiableDelayFunction):
    '''
    VDF instance over GF(p)
    '''
    
    def _setup_walk(self) :
        '''
        Prepare the isogeny walk and evaluate the basis point.

        The isogeny walk is deterministically computed from the curve
        in the trusted setup.

        OUTPUT:
        
        Triplet (E, P, data).

        - E is the image curve of the walk,
        - P is the image of the generator under the isogeny,
        - data is the evaluation key describing the walk.
        '''

        # start walking
        fP = self.P
        ek = [self.setup.E0]
        for i in range(self.delay):
            E, (fP,) = ek[-1].isogeny_forward((fP,))
            ek.append(E)
        E = ek.pop()
        ek.reverse()

        return E, fP, ek

    def evaluate(self, Q):
        'Evaluate the VDF at point Q'
        for E in self.walk_data:
            (Q,) = E.isogeny_backward(Q)
        return Q

    def verify(self, Q, fQ):
        'Verify that fQ = VDF(Q)'

        if not fQ in self.setup.E0 or fQ.is_zero() or not (self.setup.N*fQ).is_zero():
            return False

        E0 = self.setup.E0.to_gfp2()
        E1 = self.E1.to_gfp2()
        e1 = self.P.weierstrass(E0).weil_pairing(fQ.weierstrass(E0), self.setup.N)
        e2 = self.fP.weierstrass(E1).weil_pairing(Q.weierstrass(E1), self.setup.N)
        
        # # Convert to Weierstrass
        # P_ws = self.P.weierstrass()
        # fP_ws = self.fP.weierstrass()
        # Q_ws = Q.weierstrass()
        # fQ_ws = fQ.weierstrass()

        # _Z, mil11 = pairing.miller(fQ_ws, P_ws, self.setup.N, denominator=False)
        # e1 = pairing.exponentiation(self.setup, mil11[0]/mil11[1])
        # logging.debug('f_{N, fQ}(P) = %r', mil11)
        # logging.debug('e(fQ, P) = %r', e1)

        # _Z, mil22 = pairing.miller(Q_ws, fP_ws, self.setup.N, denominator=False)
        # e2 = pairing.exponentiation(self.setup, mil22[0]/mil22[1])
        # logging.debug('f_{N, Q}(fP) = %r', mil22)
        # logging.debug('e(Q, fP) = %r', e2)

        return e1 != 1 and (e1 == e2 or e1 == 1/e2)


class VDF_GFp2(VerifiableDelayFunction):
    '''
    VDF instance over GF(p²)
    '''
        
    def _setup_walk(self) :
        '''
        Prepare the isogeny walk and evaluate the basis point.

        The isogeny walk is a random walk starting from the curve
        in the trusted setup.

        OUTPUT:
        
        Triplet (E, P, data).

        - E is the image curve of the walk,
        - P is the image of the generator under the isogeny,
        - data is the evaluation key describing the walk.
        '''

        # Lift the start curve to GF(p²)
        E = self.setup.E0.to_gfp2()
        
        # start walking
        fP = self.P
        ek = [E]
        for i in range(self.delay):
            E, (fP,) = ek[-1].isogeny_forward((fP,))
            ek.append(E)
        E = ek.pop()
        ek.reverse()

        return E, fP, ek

    def evaluate(self, Q):
        'Evaluate the VDF at point Q'
        for E in self.walk_data:
            (Q,) = E.isogeny_backward(Q)
        return Q.trace(self.setup.E0)

    def verify(self, Q, fQ):
        'Verify that fQ = VDF(Q)'
        
        if not(fQ in self.setup.E0) or fQ.is_zero() or not (self.setup.N*fQ).is_zero():
            return False

        E0 = self.setup.E0.to_gfp2()
        e1 = self.P.weierstrass(E0).weil_pairing(fQ.weierstrass(E0), self.setup.N)
        e2 = self.fP.weierstrass(self.E1).weil_pairing(Q.weierstrass(self.E1), self.setup.N)**2

        # _Z, mil11 = pairing.miller(fQ_ws, P_ws, self.setup.N, denominator=False)
        # e1 = pairing.exponentiation(self.setup, mil11[0]/mil11[1])
        # logging.debug('f_{N, fQ}(P) = %r', mil11)
        # logging.debug('e(fQ, P) = %r', e1)

        # _Z, mil22 = pairing.miller(Q_ws, fP_ws, self.setup.N, denominator=True)
        # e2_squared = pairing.exponentiation(self.setup, mil22[0]/mil22[1])**2
        # logging.debug('f_{N, Q}(fP) = %r', mil22)
        # logging.debug('e(Q, fP)² = %r', e2_squared)

        return e1 != 1 and (e1 == e2 or e1 == 1/e2)

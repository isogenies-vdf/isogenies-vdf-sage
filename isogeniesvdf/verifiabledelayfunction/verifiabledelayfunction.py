# -*- coding: utf-8 -*-
import curve, point, tate
from sage.schemes.elliptic_curves.constructor import EllipticCurve
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
        self.E0_fp2 = self.setup.E0.to_gfp2()
        self.E0_fp2_ws = self.E0_fp2.weierstrass()
        self.P_ws = self.P.get_coordinates(self.E0_fp2_ws)
        self.E1_fp2 = self.E1.to_gfp2()
        self.E1_fp2_ws = self.E1_fp2.weierstrass()
        self.fP_ws = self.fP.get_coordinates(self.E1_fp2_ws)

    def __repr__(self):
        return 'Verifiable delay function with delay %d' % self.delay

    __str__ = __repr__

    def random_input(self):
        'Sample a random input point for the VDF'
        return self.E1.point_of_order(N=True, n=0, twist=False, deterministic=False)


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
            print(E)
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

        EE0 = self.E0_fp2_ws
        PP =   self.P_ws
        fQQ =  fQ.get_coordinates(EE0)
        
        EE1 = self.E1_fp2_ws
        fPP =  self.fP_ws
        QQ =  Q.get_coordinates(EE1)

        pp = tate.Tate(0)        
        e1 = pp.tate(fQQ, PP, self.setup, denominator=False)
        e2 = pp.tate(QQ, fPP, self.setup, denominator=False)
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
            E, (fP,) = ek[-1].isogeny_forward((fP,), principal=False)
            ek.append(E)
            print(E)
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
        if not(fQ in self.setup.E0.to_gfp2()) or fQ.is_zero() or not (self.setup.N*fQ).is_zero():
            return False

        pp = tate.Tate(0)
        EE0 = self.E0_fp2_ws
        PP =   self.P_ws
        fQQ =  fQ.get_coordinates(EE0)
        e1 = pp.tate(fQQ, PP, self.setup, denominator=True)
        # it should be fQQ.tate_pairing(PP, self.setup.N, 2)
        e1bis = pp.tate(2*fQQ, PP, self.setup, denominator=True)
        
        EE1 = self.E1_fp2_ws
        fPP =  self.fP_ws
        QQ =  Q.get_coordinates(EE1)
        e2 = pp.tate(fPP, QQ, self.setup, denominator=True)**2

        # on the special initial curve we could set denominator=False
        return e1 != 1 and (e1 == e2 or e1 == 1/e2)


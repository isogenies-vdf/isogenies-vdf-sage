# -*- coding: utf-8 -*- 
import pairing
from sage.rings.integer_ring import ZZ

def vdf_verif(curve, setup, Q, hat_phiQ) :
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
    
    
    #mil1 = hat_phiQ_ws._miller_(P_ws, ZZ(curve.N))
    _Z, mil11 = pairing.miller(hat_phiQ_ws, P_ws, ZZ(curve.N), denominator=False)
    e1 = pairing.exponentiation(curve, mil11[0]/mil11[1])
    #assert e1 == hat_phiQ_ws.tate_pairing(P_ws, ZZ(curve.N), 2)

    #mil2 = Q_ws._miller_(phiP_ws, ZZ(curve.N))
    _Z, mil22 = pairing.miller(Q_ws, phiP_ws, ZZ(curve.N), denominator=False)
    e2 = pairing.exponentiation(curve, mil22[0]/mil22[1])
    #assert e2 == Q_ws.tate_pairing(phiP_ws, ZZ(curve.N), 2)

    if e1 != 1 :
        if e1 == e2 :
            return True
        if e1 == 1/e2:
            return True
        # Pairing equation does not hold
        return False
    # e_Tr_hat_phiQ_P = 1
    return False

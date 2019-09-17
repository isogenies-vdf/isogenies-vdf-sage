# -*- coding: utf-8 -*- 
import pairing
from sage.rings.integer_ring import ZZ

def vdf_verif(c, setup, Q, Tr_hat_phiQ) :
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

    '''working
    print('weil')
    print(Tr_hat_phiQ_ws.weil_pairing(P_ws, ZZ(c.N)))
    print(Q_ws.weil_pairing(phiP_ws, ZZ(c.N))**2)
    print(Q_ws.weil_pairing(phiP_ws, ZZ(c.N))**(-2))

    print('tate')
    print(Tr_hat_phiQ_ws.tate_pairing(P_ws, ZZ(c.N), 2))
    print(Q_ws.tate_pairing(phiP_ws, ZZ(c.N), 2)**2)
    print(Q_ws.tate_pairing(phiP_ws, ZZ(c.N), 2)**(-2))
    '''
    
    '''not working
    print('ate')
    print(Tr_hat_phiQ_ws.ate_pairing(P_ws, ZZ(c.N), 2, ZZ(-2*c.p)))
    print(Q_ws.ate_pairing(phiP_ws, ZZ(c.N), 2, ZZ(-2*c.p))**2)
    print(Q_ws.ate_pairing(phiP_ws, ZZ(c.N), 2, ZZ(-2*c.p))**(-2))
    '''

    #mil1 = Tr_hat_phiQ_ws._miller_(P_ws, ZZ(c.N))
    _Z, mil11 = pairing.miller(Tr_hat_phiQ_ws, P_ws, ZZ(c.N), denominator=True)
    #assert mil1**((ZZ(c.p)**2-1)/ZZ(c.N)) == Tr_hat_phiQ_ws.tate_pairing(P_ws, ZZ(c.N), 2)
    e1 = pairing.exponentiation(c, mil11[0]/mil11[1])
    #assert e1 == Tr_hat_phiQ_ws.tate_pairing(P_ws, ZZ(c.N), 2)

    #mil2 = Q_ws._miller_(phiP_ws, ZZ(c.N))
    _Z, mil22 = pairing.miller(Q_ws, phiP_ws, ZZ(c.N), denominator=True)
    #assert mil2 == mil22[0]/mil22[1]
    e2_squared = pairing.exponentiation(c, mil22[0]/mil22[1])**2
    #assert e2_squared == Q_ws.tate_pairing(phiP_ws, ZZ(c.N), 2)**2

    if e1 != 1 :
        if e1 == e2_squared :
            return True
        if e1 == 1/e2_squared:
            return True
        # Pairing equation does not hold
        return False
    # e_Tr_hat_phiQ_P = 1
    return False


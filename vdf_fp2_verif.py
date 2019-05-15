# -*- coding: utf-8 -*- 
from pairing import *

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
        print 'evaluation step does not give point of the curve defined over Fp'
        return False
    
    #
    #for R in Tr_hat_phiQ :
    #    #if not(R.in_curve() and  R.x in c.Fp and R.z in c.Fp) :
    #        #print 'evaluation step does not give point of the curve defined over Fp'
    #        #return False
    #

    # this does not depend on the eval answer, can be computed before the eval
    P_ws = P.weierstrass()
    phiP_ws = phiP.weierstrass()
    Q_ws = Q.weierstrass()
    Tr_hat_phiQ_ws = Tr_hat_phiQ.weierstrass()
    
    #TODO: Change with an efficient ate pairing !
    
    phiP_miller_Q = phiP_ws._miller_(Q_ws, -2  * ZZ(c.p)-1)
    e2 = exponentiation(c, phiP_miller_Q)**2
    
    #
    #ate ?
    #Q_miller_phiP= Q_ws._miller_(phiP_ws, ZZ(c.N))
    #e2 = exponentiation(c, Q_miller_phiP)**2
    #
    
    e_Tr_hat_phiQ_P = ATE(c, Tr_hat_phiQ_ws, P_ws)
    
    #assert e_Tr_hat_phiQ_P == Tr_hat_phiQ_ws.ate_pairing(P_ws, ZZ(c.N),2, -2*ZZ(c.p))
    
    if e_Tr_hat_phiQ_P != 1 :
        if e_Tr_hat_phiQ_P == e2 :
            return True
        if e_Tr_hat_phiQ_P == 1/e2:
            return True
        print 'Pairing equation does not hold.'
        return False
    print 'e_Tr_hat_phiQ_P EQUALS 1'
    return False
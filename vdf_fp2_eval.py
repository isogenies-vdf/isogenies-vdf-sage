# -*- coding: utf-8 -*- 
from sage.all import *
import point
import curve

def vdf_eval(c, setup, Q, verbose, method):
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
        print 'Delta is not a multiple of n'
        return False
    k = ZZ(c.n//2)
    
    T = Q
    c_t = copy(c_prime)
    
    if method == 'kernel4k' :
        for R in kernelsOfBigSteps :
            R = R.change_iso_curve(c_t.a)
            [T, kernelPoint, listOfCurves] = R.isogeny_degree4k(T, k, method='withoutKernel')
            c_t = T.curve
    elif method == 'kernel4' :
        for c1 in curvesPath:
            R = point.Point(1, 1, c1).change_iso_curve(c_t.a)
            [T] = R.isogeny_degree4([T])
            c_t = T.curve

    T = T.change_iso_curve(c.a)
    #T = hatphi(Q)
    
    #assert montgomery_ladder(N, T, a)[1] == 0
    
    #Trace trick
    frob_T = point.Point(T.x**c.p, T.z**c.p, T.curve)
    
    #assert montgomery_ladder(N, frob_T, a)[1] == 0
        
    #not efficient
    fQ_ws = frob_T.weierstrass()
    Q_ws = T.weierstrass()
    
    R1 = Q_ws + fQ_ws
    R2 = Q_ws - fQ_ws
    
    #the (+/-) point to return is the one defined over Fp :-)
    return c.getPointFromWeierstrass(R1) if R1[0] in c.Fp and R1[1] in c.Fp else c.getPointFromWeierstrass(R2)
    #return [c.getPointFromWeierstrass(R1), c.getPointFromWeierstrass(R2)]

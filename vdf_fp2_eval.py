# -*- coding: utf-8 -*- 
from point import Point
import curve
from sage.rings.integer_ring import ZZ

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
            [T, kernelPoint, listOfCurves] = R.isogeny_degree4k_strategy(T, k, method='withoutKernel')
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

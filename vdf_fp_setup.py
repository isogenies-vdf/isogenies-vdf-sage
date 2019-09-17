# -*- coding: utf-8 -*- 
from point import Point
from sage.rings.integer_ring import ZZ
from copy import copy
from sage.misc.misc import cputime

def isogeny_walk(curve, P, verbose, method) :
    '''
    INPUT:
    * curve the elliptic curve
    * P a point of the initial elliptic curve
    * verbose for the comments
    * method for the method of storing the isogeny (kernel4, kernel4k or withoutKernel)
    OUTPUT:
    * curve_prime the final elliptic curve of the protocol
    * final_hat_phi the list of the dual of 4 (or 4^k)-isogenies kernels for the inverse walk, with the curves associated
    * phiP the image point of P by the isogeny
    REMARK:
    * Do not walk to the j=0,1728 curve at first step! The formulas of [4]-isogenies for these curves do not work for the moment.
    '''
    
    # Delta = (n-2) * q
    # q isogenies of degree 2^(n-2)=4^(n//2 - 1)    
    q = curve.Delta // (curve.n - 2)
    if curve.Delta % (2*curve.n-4) != 0 :
        raise RuntimeError('Delta is not a multiple of 2*(n-2)')
    
    #we use a n-1 order point to define a n-2 step, in order to stay on the cratere
    k = ZZ(curve.n-1)//2
    
    ev_P = P
    curvesPath = []
    kernelsOfBigSteps = []
    
    curve_prime = copy(curve)
    
    for i in range(q) :
        #P4k defines the kernel of the isogeny
        #Attention, we need to choose the right direction ! With the twist, we go to j=0 or 1728 curve, and there is a problem with the formulas of 4-isogeny...
        P4k = curve_prime.power_of_2_order_random_point(2*k, 1, False)
        #Warning ! We do a 4**(k-1) isogeny !
        [ev_P, kernelDual, listOfCurves] = P4k.isogeny_degree4k_strategy(ev_P, k, method, stop=1)
        if not(kernelDual is None) :
            kernelsOfBigSteps += [kernelDual]
        if not(listOfCurves is None) :
            curvesPath += listOfCurves
        curve_prime = ev_P.curve
    phiP = ev_P
    curvesPath = curvesPath[::-1]
    kernelsOfBigSteps = kernelsOfBigSteps[::-1]
    return [curve_prime, curvesPath, kernelsOfBigSteps, phiP]

def vdf_setup(curve, verbose, method) :
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

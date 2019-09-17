# -*- coding: utf-8 -*- 
from sage.all import *
from point import Point

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
    * phiP the image point by the isogeny
    REMARK:
    * We do not walk to the j=0,1728 curve at first step!
    '''

    # Delta = n * q
    # q isogenies of degree 2^n=4^(n//2)
    q = curve.Delta // curve.n
    if curve.Delta % curve.n != 0 :
        raise Runtime('Delta is not a multiple of n')

    #we do 4-isogenies and not 2-isogenies
    k = ZZ(curve.n//2)

    ev_P = P
    curvesPath = []
    kernelsOfBigSteps = []

    curve_prime = copy(curve)
    first = True
    for i in range(q) :
        #P4k defines the kernel of the isogeny
        if first : # we check the first step does not go to j=0 or 1728 curve
            j = 0
            while j == 0 or j == 1728 :
                P4k = curve_prime.power_of_2_order_random_point(2*k, 2)
                P4 = P4k.get_P4(k)
                xP4 = P4.normalize().x
                while xP4 == 1 or xP4 == -1 :
                    P4k = curve_prime.power_of_2_order_random_point(2*k, 2)
                    P4 = P4k.get_P4(k)
                    xP4 = P4.normalize().x
                #assert isOrder4k(k, P4k, aprime)
                curve_onestep = P4.isogeny_degree4([ev_P])[0].curve
                j = curve_onestep.weierstrass().j_invariant()
            first = False
        else : # we ckeck that the step does not backtrack
            P4k = curve_prime.power_of_2_order_random_point(2*k, 2)
            P4 = P4k.get_P4(k)
            xP4 = P4.normalize().x
            while xP4 == 1 or xP4 == -1 :
                P4k = curve_prime.power_of_2_order_random_point(2*k, 2)
                P4 = P4k.get_P4(k)
                xP4 = P4.normalize().x
        [ev_P, kernelDual, listOfCurves] = P4k.isogeny_degree4k_strategy(ev_P, k, method)
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


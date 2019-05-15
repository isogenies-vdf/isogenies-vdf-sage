# -*- coding: utf-8 -*- 
from sage.all import *
import curve

def double_line(P, Q, a) :
    # P is in jacobian coordinates and Q in affine coordinates
    # returns the line tangent in P evaluated in Q, as a fraction, and the point 2*P in jacobian coordinates
    
    X, Y, Z = P[0], P[1], P[2]
    xQ, yQ = Q[0], Q[1]
    
    t1 = Y**2
    
    t2 = X * t1
    t2 = t2 + t2
    t2 = t2 + t2
    # t2 = 4 * X * Y**2
    
    t3 = t1**2
    t3 = t3 + t3
    t3 = t3 + t3
    t3 = t3 + t3
    # t3 = 8 * Y**4
    
    if len(P) == 4 :
        t4 = P[3]
    else : #for the non-new curves... 
        t4 = Z**2
        
    if a == -3 :
        t5 = (X - t4) * (X + t4)
        t5 = t5 + t5 + t5
    else :
        t5 = X**2
        t5 = t5 + t5 + t5
        t6 = t4**2
        t5 = t5 + a * t6
    # t5 = 3 * X**2 + a * Z**4
    
    XX = t5**2 - t2 - t2
    # XX = (3*X**2 + a*Z**4)**2 - 8*X*Y**2
    YY = t5 * (t2 - XX) - t3
    # YY = (3*X**2 + a*Z**4)*(4*X*Y**2 - XX) - 8*Y**4
    ZZ = Y * Z
    ZZ = ZZ + ZZ
    # ZZ = 2*Y*Z
    TT = ZZ**2
    
    # 6Spk + 3Mpk if a petit != -3
    # 4Spk + 4Mpk if a == -3
    
    
    mu_num = t5
    # mu_num = 3 * X**2 + a * Z**4
    
    #
    l_den = ZZ * t4
    #
    
    # l_den = 2 * Y * Z**3
    l_num = yQ * l_den - t1 - t1 - t5 * (t4 * xQ - X)
    # l_num = yQ*l_den - 2 * Y**2 - mu_num * (Z**2* xQ - X)
    
    # 2 Mpk + 2k Mp
    
    
    # Total:
    #  if a == -3 : 4 Spk + 6 Mpk + 2k Mp
    #  if a != -3 : 6 Spk + 5 Mpk + 2k Mp
    
    return [[l_num, l_den], [XX, YY, ZZ, TT]]

def add_line(R, P, Q) :
    # R is in jacobian coordinates and P, Q in affine coordinates
    # returns the line P-R evaluated in Q, as a fraction, and the point P+R in jacobian coordinates
    
    X1, Y1, Z1 = R[0], R[1], R[2]
    X, Y = P[0], P[1]
    xQ, yQ = Q[0], Q[1]
    
    if len(R) == 4 :
        t1 = R[3]
    else : #useful for non-new curves...
        t1 = Z1**2
    t2 = Z1 * t1        # Z1**3
    t3 = X * t1         # XZ1**2
    t4 = Y * t2         # YZ1**3
    t5 = t3 - X1        # XZ1**2 - X1
    t6 = t4 - Y1        # YZ1**3 - Y1
    t7 = t5**2           # (XZ1**2-X1)**2
    t8 = t5 * t7        # (XZ1**2-X1)**3
    t9 = X1 * t7        # (XZ1**2-X1)**2Z1**2
    
    X3 = t6**2 - (t8 + t9 + t9)
    Y3 = t6 * (t9 - X3) - Y1 * t8
    Z3 = Z1 * t5
    T3 = Z3**2
    # 3 Spk + 8 Mpk
    

    """
    X3 = (Y*Z1**3 - Y1)**2 - (X*Z1**2 - X1)**2 * (X1 + X*Z1**2)
    Y3 = (Y*Z1**3 - Y1) * ((X*Z1**2 - X1)**2*X1 - X3) - Y1 * (X*Z1**2 - X1)**3
    Z3 = (X*Z1**2 - X1) * Z1
    T3 = Z3**2
    """
    
    mu_num = t6
    # mu_num = (Y*Z1**3 - Y1)
    
    mu_den = Z3
        
    # mu_den = Z1 * (X * Z1**2 - X1)
    l_num = mu_den * (yQ - Y) - mu_num * (xQ - X)
    # l_num =  Z1 * (X * Z1**2 - X1) * (yQ - Y) - (Y*Z1**3 - Y1) * (xQ - X)
    
    l_den = mu_den
    # l_den = Z1 * (X * Z1**2 - X1)
    
    # 2 Mpk
    
    # Total: 10 Mpk + 3 Spk
    return [[l_num, l_den], [X3, Y3, Z3, T3]]

def miller(P, Q, n, a, denominator=False) :
    # return the miller loop f_{P, n}(Q) for an even embedded degree curve
    if Q.is_zero():
        raise ValueError("Q must be nonzero.")
    n = ZZ(n)
    if n.is_zero():
        raise ValueError("n must be nonzero.")
    n_is_negative = False
    if n < 0:
        n = n.abs()
        n_is_negative = True
        
    t_num, t_den = 1, 1
    V = P
    S = 2*V # V=P is in affine coordinates
    nbin = n.bits()
    i = n.nbits() - 2
    while i > -1:
        [[ell_num, ell_den], S] = double_line(V, Q, a)
        t_num =(t_num**2)*ell_num
        if denominator :
            t_den =(t_den**2)*ell_den
        
        V = S
        if nbin[i] == 1:
            [[ell_num, ell_den], S] = add_line(V, P, Q)
            t_num = t_num*ell_num
            if denominator :
                t_den = t_den*ell_den
            
            V = S
        i = i-1
    
    if not(denominator) :
        t_den = 1
    
    if n_is_negative :
        t_num, t_den = t_den, t_num
        S[1] = -S[1]
    return [[t_num,t_den], S]

def exponentiation(curve, x) :
    """
    INPUT: 
    * x an element of Fp2
    OUTPUT:
    * x**((p**2-1)//N) in our special case where (p**2 - 1)//N == f* 2**(n+1) * (2**(n-1) * f * N - 1)
    """
    x1 = (x**curve.f)
    for i in range(curve.n+1) :
        x1 = x1**2
    x2 = x1
    x3 = x2**(curve.f*curve.N)
    for i in range(curve.n-1) :
        x3 = x3**2
    return x3/x2

def miller_ATE(curve, P, Q, denominator=False) :
    t = -2 * curve.p
    a = curve.weierstrass().a4()
    m1, L = miller(Q, P, t - 1, a, denominator)
    return [m1[0],m1[1]] 

def ATE(curve, P, Q, denominator=False) :
    t = -2 * curve.p
    a = curve.weierstrass().a4()
    m1, L = miller(Q, P, t - 1, a, denominator)
    #maybe an inversion if t-1 < 0
    m2 = m1[0]/m1[1]
    return exponentiation(curve, m2)

def TATE(curve, P, Q, denominator=False) :
    #if the curve is defined over Fp, we do not compute the denominators
    a = curve.weierstrass().a4()
    m1, L = miller(P, Q, curve.N, a, denominator)
    m2 = m1[0]/m1[1]
    return exponentiation(curve, m2)
    #if the curve is defined over Fp2, we need to compute denominators
    #my double_line and add_line do not work for the moment
    #I use the Sage ones...
    
    #R = (curve.p+1)//curve.N * Q.parent().random_point()
    #while R == 0 : 
    #    R = (curve.p+1)//curve.N * Q.parent().random_point()
    #assert P.weil_pairing(Q+R, ZZ(curve.N)) !=  P.weil_pairing(R, ZZ(curve.N))
    #m1, L = miller(P, Q+R, curve.N, a, denominator)
    #m11, L = miller(P, R, curve.N, a)
    #m2 = m1[0] * m11[1]/(m1[1]* m11[0])


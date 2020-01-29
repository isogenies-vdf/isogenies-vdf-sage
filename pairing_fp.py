#!/usr/bin/env sage
# -*- coding: utf-8 -*-

def double_line_J_montgomery(S, P, a2):
    X1, Y1, Z1 = S[0], S[1], S[2]
    xP, yP = P[0], P[1]
    
    # doubling formulas
    t1 = Y1**2          # Y1²
    t2 = 4*X1*t1        # 4 X1 Y1²
    t3 = 8*t1**2        # 8 Y1^4
    t4 = Z1**2          # Z1²
    t5 = 3*X1**2 + 2*a2*X1*t4 + t4**2 #numerator after projectivation
    Z3 = 2*Y1*Z1
    X3 = t5**2 - 2*t2 - a2*Z3**2
    Y3 = t5*(t2 - X3) - t3

    #ss = 2*E.point([X1/Z1**2, Y1/Z1**3, 1])
    #assert ss == E.point([X3/Z3**2, Y3/Z3**3, 1])

    l_den = Z3*t4       # 2 Y1 Z1³
    l_num = l_den*yP - 2*t1 + t5 * (X1 - t4*xP)

    return [[l_num, l_den], [X3,Y3,Z3]]

def add_line_J_montgomery(S, P, Q, a2):
    X1, Y1, Z1 = S[0], S[1], S[2]
    X, Y = P[0], P[1]
    assert P[2] == 1
    xQ, yQ = Q[0], Q[1]
    assert Q[2] == 1
    
    # addition formulas
    t1 = Z1**2      # Z1²
    t2 = Z1*t1      # Z1³
    t3 = X*t1       # Z1² X
    t4 = Y*t2       # Y Z1³
    t5 = t3 - X1    # Z1² X - X1
    t6 = t4 - Y1    # Y Z1³ - Y1
    t7 = t5**2      # (Z1² X - X1)²
    t8 = t5*t7      # (Z1² X - X1)³
    t9 = X1 * t7    # X1*(Z1² X - X1)²
    X3 = t6**2 - (t8+2*t9+t7*t1*a2)  # (Z1²Y-Y1)² - (Z1²X-X1)²(Z1²X+Z1²a2+X1)
    Y3 = t6*(t9 - X3) - Y1*t8
    Z3 = Z1*t5      # Z1*(Z1² X - X1)

    #assert Z1 != 0 or E.point([X1/Z1**2, Y1/Z1**3, 1]) + P == E.point([X3/Z3**2, Y3/Z3**3, 1])

    # line in the general case
    #(l_den = Z3)
    l_num = Z3 * (yQ - Y) - t6 * (xQ - X)
    
    return [[l_num, Z3], [X3, Y3, Z3]]

def miller(P, Q, n, a2, denominator=True) :
    # return the miller loop f_{P, n}(Q)
    if Q.is_zero():
        raise ValueError("Q must be nonzero.")
    if n.is_zero():
        raise ValueError("n must be nonzero.")
    n_is_negative = False
    if n < 0:
        n = n.abs()
        n_is_negative = True
        
    t_num, t_den = 1, 1
    V = (P[0], P[1], 1)
    assert P[2] == 1
    assert Q[2] == 1

    S = 2*P # P is in affine coordinates so it works
    nbin = n.bits()
    i = n.nbits() -2
    while i > -1:
        [[ell_num, ell_den], S] = double_line_J_montgomery(V, Q, a2)
        t_num = ell_num * t_num**2
        if denominator :
            t_den = ell_den * t_den**2
        V = S
        if nbin[i] == 1:
            [[ell_num, ell_den], S] = add_line_J_montgomery(V, P, Q, a2)
            t_num = t_num*ell_num
            if denominator :
                t_den = t_den*ell_den
            V = S
        i -= 1
    if not(denominator) :
        t_den = 1

    if n_is_negative :
        t_num, t_den = t_den, t_num
        S[1] = -S[1]
    return [[t_num,t_den], S]

def exponentiation(setup, x) :
    """
    INPUT: 
    * setup the setup for the vdf
    * x an element of Fp2
    OUTPUT:
    * x**((p**2-1)//N) in our special case where (p**2 - 1)//N == f* 2**(n+1) * (2**(n-1) * f * N - 1)
    """
    x1 = (x**setup.f)
    for i in range(setup.n+1) :
        x1 = x1**2
    x2 = x1
    x3 = x2**(setup.f*setup.N)
    for i in range(setup.n-1) :
        x3 = x3**2
    return x3/x2

def tate(setup, P, Q, denominator=False) :
    #if the curve is defined over Fp, we do not compute the denominators
    a2 = P.curve().a2()
    m1, L = miller(P, Q, setup.N, a2, denominator)
    m2 = m1[0]/m1[1]
    return exponentiation(setup, m2)
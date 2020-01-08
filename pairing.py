# -*- coding: utf-8 -*- 

def eval_line(R, P, Q) :
    '''
    INPUT: 
    - R a point of an EllipticCurve object
    - P a point of an EllipticCurve object
    - Q a point of an EllipticCurve object
    OUPUT:
    - The line through P and R evaluated at Q, given with the numerator and
      the denominator.
    - The point P+R
    '''
    if Q.is_zero():
        raise ValueError("Q must be nonzero.")

    if P.is_zero() or R.is_zero():
        if P == R:
            return [P, [P.curve().base_field().one(), 1]]
        if P.is_zero():
            return [R, [Q[0] - R[0], 1]]
        if R.is_zero():
            return [P, [Q[0] - P[0], 1]]
    elif P != R:
        if P[0] == R[0]:
            return [P.curve()(0), [Q[0] - P[0], 1]]
        else:
            lnum, lden = (R[1] - P[1]), (R[0] - P[0])
            t1 = lden**2
            t2 = t1*lden
            xPplusR = (lnum**2) - (lden**2)*(P.curve().a2() + P[0] + R[0])
            yPplusR = R[1]*(t2) + lnum*(xPplusR - (t1)*R[0])
            zPplusR = t2
            PplusR = P.curve()(lden*xPplusR, -yPplusR, zPplusR)
            return [PplusR, [lden * (Q[1] - P[1])  - lnum * (Q[0] - P[0]), lden]]
    else:
        a1, a2, a3, a4, a6 = P.curve().a_invariants()
        numerator = (3*P[0]**2 + 2*a2*P[0] + a4 - a1*P[1])
        denominator = (2*P[1] + a1*P[0] + a3)
        if denominator == 0:
            return [P.curve()(0), [Q[0] - P[0], 1]] #except in characteristic 2 ?
        else:
            #l = numerator/denominator
            t1 = denominator**2
            t2 = t1 * denominator
            x2P = numerator**2 - (P.curve().a2() + 2* P[0]) * t1
            y2P = P[1]*t2 + numerator*(x2P - t1 * P[0])
            z2P = t2
            twoP = P.curve()(denominator*x2P, -y2P, z2P)
            return [twoP, [denominator * (Q[1] - P[1]) - numerator * (Q[0] - P[0]), denominator]]

def miller(P, Q, n, denominator=False) :
    '''
    INPUT:
    - P a point of an EllipticCurve object
    - Q a point of an EllipticCurve object
    - n an integer (the Miller loop)
    OUPUT:
    - the point [n]P
    - the Miller loop f_{P, n}(Q) in the context of our curves, given with numerator and denominator.
    '''
    if Q.is_zero():
        raise ValueError("Q must be nonzero.")
    if n == 0:
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
        [S, [ell_num, ell_den]] = eval_line(V, V, Q)
        t_num = (t_num**2)*ell_num
        if denominator :
            [R, [vee_num, vee_den]] = eval_line(S, -S, Q)
            t_den = (t_den**2)*ell_den*vee_num
            t_num *= vee_den
        V = S
        if nbin[i] == 1:
            [S, [ell_num, ell_den]] = eval_line(V, P, Q)
            t_num = t_num*ell_num
            if denominator :
                [R, [vee_num, vee_den]] = eval_line(S, -S, Q)
                t_den *= ell_den*vee_num
                t_num *= vee_den
            V = S
        i = i-1
    if not(denominator) :
        t_den = 1
    if n_is_negative :
        t_num, t_den = t_den, t_num
        S = -S
    return [S, [t_num,t_den]]


def exponentiation(setup, x):
    '''
    INPUT: 
    * x an element of Fp2
    OUTPUT:
    * x**((p**2-1)//N) in our special case where (p**2 - 1)//N == f* 2**(n+1) * (2**(n-1) * f * N - 1)
    '''
    x1 = (x**setup.f)
    for i in range(setup.n+1) :
        x1 = x1**2
    x2 = x1
    x3 = x2**(setup.f*setup.N)
    for i in range(setup.n-1) :
        x3 = x3**2
    return x3/x2

def tate(P, Q, E, denominator=True):
    '''
    INPUT:
    - P a point of an EllipticCurve object
    - Q a point of an EllipticCurve object
    - E an elliptic curve object
    OUPUT:
    - the Tate pairing T(P, Q)
    '''
    mill_nd = miller(P, Q, E.setup.N, denominator=denominator)[1]
    return exponentiation(E.setup, mill_nd[0] / mill_nd[1])

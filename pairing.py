# -*- coding: utf-8 -*- 
from copy import copy

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
        if P[0]*R[2] == R[0] * P[2] and P[1]*R[2] == R[1] * P[2]:
            return [P, [P.curve().base_field().one(), 1]]
        if P.is_zero():
            return [R, [R[2]* Q[0] - Q[2] * R[0], Q[2] * R[2]]]
        if R.is_zero():
            return [P, [P[2] * Q[0] - Q[2] * P[0], Q[2] * P[2]]]
    elif P != R:
        if P[0]*R[2] == P[2] * R[0]:
            return [P.curve().point([0,1,0], check=False), [P[2] * Q[0] - Q[2] * P[0], Q[2] * P[2]]]
        else:
            n,d = (P[2] * R[1] - R[2] * P[1]), (P[2] * R[0] - R[2] * P[0])
            '''
            Formulas with normalization :
            X = (n/d)² - (a+xp/zp + xr/zr)
            Y = -(n/d * xpplusr + (d*yr - n*xr)/(d*zr))
            Z = 1
            and line
            yq/zq - n/d * xq/zq - (d*yr-n*xr)/(d*zr)

            which gives:

            X = n² * zp * zr - d² * (zp*zr*a + zr*xp + zp*xr)
            Y = -(n*xpplusr/d + d*zp*(d*yr - n*xr)) # TODO one inversion here
            Z = d²*zp*zr
            line_n = d*zr*yq - n*zr*xq - zq * (d*yr - n*xr)
            line_d = d*zr*zq
            '''
            t1 = d**2
            t2 = d*R[2]
            t3 = n*R[0]
            X = n**2 * P[2] * R[2] - t1 * (P[2] * R[2] * P.curve().a2() + R[2]*P[0] + P[2]*R[0])
            Y = -(n*X/d + d*P[2]*(d*R[1] - t3)) # TODO one inversion here
            Z = t1 * P[2] * R[2]
            PplusR = P.curve().point([X,Y,Z], check=False)
            line_n = t2 * Q[1] - n * R[2]*Q[0] - Q[2]*(d*R[1] - t3)
            line_d = t2*Q[2]
            return [PplusR, [line_n, line_d]]
    else:
        a1, a2, a3, a4, a6 = P.curve().a_invariants()
        t1 = P[0] * P[2]
        t2 = P[1] * P[2]
        t3 = P[2]**2
        n = (3*P[0]**2 + 2*a2*t1 + a4*t3 - a1*t2)
        d = 2*t2 + a1*t1 + a3*t3
        if d == 0:
            return [P.curve().point([0,1,0], check=False), [P[2] * Q[0] - P[0]*Q[2], Q[2]* P[2]]] #except in characteristic 2 ?
        else:
            '''
            Formulas with normalization :
            X = (n/d)² - (a+2*xp/zp)
            Y = -(n/d * xpplusr + (d*yp - n*xp)/(d*zp))
            Z = 1
            and line
            yq/zq - n/d * xq/zq - (d*yp-n*xp)/(d*zp)

            which gives:

            X = n² * zp - d² * (zp*a + 2*xp)
            Y = -(n*X/d + d*(d*yp - n*xp)) # TODO one inversion here
            Z = d²*zp
            line_n = d*zp*yq - n*zp*xq - zq * (d*yp - n*xp)
            line_d = d*zp*zq
            '''
            
            t1 = d**2
            t2 = d*P[2]
            t3 = n*P[0]
            X = n**2 * P[2] - t1 * (P[2] * P.curve().a2() + 2*P[0])
            Y = -(n*X/d + d*(d*P[1] - t3)) # TODO one inveresion here
            Z = t1*P[2]
            twoP = P.curve().point([X,Y,Z], check=False)
            line_n = t2 * Q[1] - n * P[2] * Q[0] - Q[2] *(d*P[1] - t3)
            line_d = t2 * Q[2]
            return [twoP, [line_n, line_d]]

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

    # S = 2*V computed in projective coordinates formula (not done in sage)
    a1, a2, a3, a4, a6 = V.curve().a_invariants()
    ell_n = (3*V[0]**2 + 2*a2*V[0]*V[2] + a4*V[2]**2 - a1*V[1]*V[2])
    ell_d = 2*V[1]*V[2] + a1*V[0]*V[2] + a3*V[2]**2
    X = ell_n**2 * V[2] - ell_d**2 * (V[2] * P.curve().a2() + 2*V[0])
    Y = -(ell_n*X/ell_d + ell_d*(ell_d*V[1] - ell_n*V[0])) # TODO one inveresion here
    Z = ell_d**2*V[2]
    S = V.curve().point([X,Y,Z], check=False)
    
    nbin = n.bits()
    i = n.nbits() - 2
    while i > -1:
        [S, [ell_num, ell_den]] = eval_line(V, V, Q)
        t_num = (t_num**2)*ell_num
        if denominator :
            minusS = S.curve().point([S[0], -S[1], S[2]], check=False)
            [R, [vee_num, vee_den]] = eval_line(S, minusS, Q)
            t_den = (t_den**2)*ell_den*vee_num
            t_num *= vee_den
        V = S
        if nbin[i] == 1:
            [S, [ell_num, ell_den]] = eval_line(V, P, Q)
            t_num = t_num*ell_num
            if denominator :
                minusS = S.curve().point([S[0], -S[1], S[2]], check=False)
                [R, [vee_num, vee_den]] = eval_line(S, minusS, Q)
                t_den *= ell_den*vee_num
                t_num *= vee_den
            V = S
        i = i-1
    if not(denominator) :
        t_den = 1
    if n_is_negative :
        t_num, t_den = t_den, t_num
        S[1] *= -1
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

#!/usr/bin/env sage
# -*- coding: utf-8 -*-

import setup, curve, point, vdf, pairing
from sage.schemes.elliptic_curves.constructor import EllipticCurve

def test_isog(E):
    '''
    Test that composing isogenies forth and back equals mul by 4.
    '''
    Ps = [E.elligator(2), E.elligator(10), E.elligator(111)]
    E1, Qs = E.isogeny_forward(Ps)
    for A, B in zip(E.isogeny_backward(*Qs), Ps):
        assert A == 2*B, 'Curve %r, %r != %r' % (E, A.normalize(), B.normalize())
    return E1

def test_vdf(setup, delay, vdfclass=vdf.VDF_GFp, reps=10):
    vdf = vdfclass(setup, delay)
    for i in range(reps):
        Q = vdf.random_input()
        fQ = vdf.evaluate(Q)
        assert vdf.verify(Q, fQ)

def test_double_line_jac(E, reps=10):
    EE = E.weierstrass()
    for i in range(reps):
        S = EE.random_point()
        while S.is_zero():
            S = EE.random_point()
        P = EE.random_point()
        while P.is_zero():
            P = EE.random_point()
        u = E.field.random_element()
        S_j = (S[0]*u**2, S[1]*u**3, u)
        [l, dblS] = pairing.double_line_jac(S_j, P, EE.a2())
        assert (dblS[2] == 0 and S[1] == 0) or 2*S == EE.point([dblS[0]/dblS[2]**2, dblS[1]/dblS[2]**3])
        assert S._line_(S, P) == l[0]/l[1]

def test_add_line_jac(E, reps=10):
    EE = E.weierstrass()
    for i in range(reps):
        S = EE.random_point()
        while S.is_zero():
            S = EE.random_point()
        P = EE.random_point()
        while P.is_zero():
            P = EE.random_point()
        Q = EE.random_point()
        while Q.is_zero():
            Q = EE.random_point()
        u = E.field.random_element()
        S_j = (S[0]*u**2, S[1]*u**3, u)
        [l, SplusP] = pairing.double_line_jac(S_j, P, Q, EE.a2())
        assert S+P == EE.point([SplusP[0]/SplusP[2]**2, SplusP[1]/SplusP[2]**3])
        assert S._line_(P, Q) == l[0]/l[1]

def test_miller(E, reps=10):
    E_Fp2 = E.to_gfp2()
    EE = EllipticCurve(E_Fp2.field, [0,E_Fp2.A,0,1,0])
    i = 0
    p = E_Fp2.field.characteristic()
    while i<reps:
        P =  E.point_of_order(N=True, n=0, twist=False, deterministic=False)
        PP = P.get_coordinates(EE)
        PP_j = (PP[0]*PP[2], PP[1]*PP[2]**2, PP[2])
        Q = E_Fp2.point_of_order(N=True, n=0, twist=False, deterministic=False)
        while Q.x/Q.z in E.field:
            Q = E_Fp2.point_of_order(N=True, n=0, twist=False, deterministic=False)
        QQ = Q.get_coordinates(EE)
        assert (E.setup.N * P).is_zero() and not(P.is_zero())
        assert (E.setup.N * Q).is_zero() and not(Q.is_zero())
        [l1, S1] = pairing.miller(PP_j, QQ, E.setup.N, EE.a2(), denominator=True)
        l2 = PP._miller_(QQ, E.setup.N)
        assert l2 == l1[0]/l1[1]
        assert S1[2] == 0
        i+=1

def test_sqrt(curve, reps=10):
    fp2 = curve.to_gfp2().field
    list_squares = [(fp2.random_element())**2 for i in range(reps)]
    for i in range(reps):
        sqrt_u = curve._sqrt(list_squares[i])
        assert sqrt_u**2 == list_squares[i]

def test_exponentiation(setup, reps=10):
    F = setup.E0.field
    for i in range(reps):
        x = F.random_element()
        assert pairing.exponentiation(setup, x) == x**((F.characteristic()**2 - 1)//setup.N)

def test_tate(E, reps=10):
    E_Fp2 = E.to_gfp2()
    EE = E_Fp2.weierstrass()
    i = 0
    p = E_Fp2.field.characteristic()
    while i<reps:
        P =  E.point_of_order(N=True, n=0, twist=False, deterministic=False)
        PP = P.get_coordinates(EE)
        Q = E_Fp2.point_of_order(N=True, n=0, twist=False, deterministic=False)
        while Q.x/Q.z in E.field:
            Q = E_Fp2.point_of_order(N=True, n=0, twist=False, deterministic=False)
        QQ = Q.get_coordinates(EE)
        ePQ = pairing.tate(E.setup, PP, QQ)
        assert ePQ == PP.tate_pairing(QQ, E.setup.N, k=2)
        i+=1

    # Check the pairing computation without denominators
    if E.field.is_prime_field() :
        i = 0
        while i<reps:
            P =  E.point_of_order(N=True, n=0, twist=False, deterministic=False)
            PP = P.get_coordinates(EE)
            Q = E.point_of_order(N=True, n=0, twist=True, deterministic=False)
            #while Q.x/Q.z in E.field:
            #    Q = E_Fp2.point_of_order(N=True, n=0, twist=False, deterministic=False)
            QQ = Q.get_coordinates(EE)
            ePQ = pairing.tate(E.setup, PP, QQ, denominator=False)
            assert ePQ == PP.tate_pairing(QQ, E.setup.N, k=2)
            i+=1

if __name__ == '__main__':
    for key, s in setup.SETUPS.items():
        print('\n==== %s ====' % key)
        print('Testing isogeny formula')
        E1 = test_isog(s.E0)
        test_isog(E1)
        print('Testing GF(p) VDF')
        test_vdf(s, 100)
        print('Testing GF(p²) VDF')
        test_vdf(s, 10, vdf.VDF_GFp2)
        print('Testing double line evaluation in jacobian coordinates')
        test_double_line_jac(s.E0)
        test_double_line_jac(E1)
        print('Testing add line evaluation in jacobian coordinates')
        test_double_line_jac(s.E0)
        test_double_line_jac(E1)
        print('Testing miller loop')
        test_miller(s.E0)
        test_miller(E1)
        print('Testing fast exponentiation')
        test_exponentiation(s)
        print('Testing Tate pairing')
        test_tate(s.E0)
        test_tate(E1)
        print('Testing efficient square root')
        test_sqrt(s.E0)
        test_sqrt(E1)


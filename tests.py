#!/usr/bin/env sage
# -*- coding: utf-8 -*-

import setup, curve, point, vdf, pairing
from sage.schemes.elliptic_curves.constructor import EllipticCurve

import sage.structure.proof.all as prf
prf.arithmetic(False)

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

def test_eval_line(E, reps=10):
    it = E.iterelligator(twist=False, deterministic=False)
    EE = EllipticCurve(E.field, [0,E.A,0,1,0])
    i = 0
    while i<reps:
        P = next(it)
        PP = P.get_coordinates(EE)
        Q = next(it)
        QQ = Q.get_coordinates(EE)
        Sum,line = pairing.eval_line(PP, QQ, QQ)
        assert line[0] == 0
        assert Sum == PP + QQ
        Sum2,line2 = pairing.eval_line(PP, PP, PP)
        assert line2[0] == 0
        assert Sum2 == 2*PP
        i+=1

def test_miller(E, reps=10):
    E_Fp2 = E.to_gfp2()
    EE = EllipticCurve(E_Fp2.field, [0,E_Fp2.A,0,1,0])
    i = 0
    p = E_Fp2.field.characteristic()
    while i<reps:
        P =  E.point_of_order(N=True, n=0, twist=False, deterministic=False)
        PP = P.get_coordinates(EE)
        Q = E_Fp2.point_of_order(N=True, n=0, twist=False, deterministic=False)
        while Q.x/Q.z in E.field:
            Q = E_Fp2.point_of_order(N=True, n=0, twist=False, deterministic=False)
        QQ = Q.get_coordinates(EE)
        assert (E.setup.N * P).is_zero() and not(P.is_zero())
        assert (E.setup.N * Q).is_zero() and not(Q.is_zero())
        S1, [n1,d1] = pairing.miller(PP, QQ, s.N, denominator=True)
        S2, [n2,d2] = pairing.miller(2*PP, QQ, s.N, denominator=True)
        S3, [n3,d3] = pairing.miller(PP, 2*QQ, s.N, denominator=True)
        e2PQ = (n2/d2)**((p**2-1)//s.N)
        eP2Q = (n3/d3)**((p**2-1)//s.N)
        ePQ2 = (n1/d1)**((p**2-1)//s.N*2)
        assert e2PQ == eP2Q and eP2Q == ePQ2
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
    EE = EllipticCurve(E_Fp2.field, [0,E_Fp2.A,0,1,0])
    i = 0
    p = E_Fp2.field.characteristic()
    while i<reps:
        P =  E.point_of_order(N=True, n=0, twist=False, deterministic=False)
        PP = P.get_coordinates(EE)
        Q = E_Fp2.point_of_order(N=True, n=0, twist=False, deterministic=False)
        while Q.x/Q.z in E.field:
            Q = E_Fp2.point_of_order(N=True, n=0, twist=False, deterministic=False)
        QQ = Q.get_coordinates(EE)
        assert pairing.tate(PP, 2*QQ, E, denominator=True) == pairing.tate(2*PP, QQ, E, denominator=True)
        assert pairing.tate(PP, QQ, E, denominator=True)**2 == pairing.tate(2*PP, QQ, E, denominator=True)
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
            assert pairing.tate(PP, 2*QQ, E, denominator=True) == pairing.tate(2*PP, QQ, E, denominator=False)
            assert pairing.tate(PP, QQ, E, denominator=True)**2 == pairing.tate(2*PP, QQ, E, denominator=False)
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
        print('Testing line evaluation')
        test_eval_line(s.E0)
        test_eval_line(E1)
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

#!/usr/bin/env sage
# -*- coding: utf-8 -*-

import setup, curve, point, vdf, pairing
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.misc.misc import cputime

import sage.structure.proof.all as prf
prf.arithmetic(False)

def bench_sqrt(curve, reps=10):
    '''
    Compute a list of reps squares and benchmark the sqrt() and our _sqrt
    functions.
    '''
    fp2 = curve.to_gfp2().field
    list_squares = [(fp2.random_element())**2 for i in range(reps)]
    t1 = cputime()
    for i in range(reps):
        sqrt_u = list_squares[i].sqrt()
    t1 = cputime(t1)
    print('Sage version:\t%s s\t(%s reps)' % (t1, reps))
    t2 = cputime()
    for i in range(reps):
        sqrt_u = curve._sqrt(list_squares[i])
    t2 = cputime(t2)
    print('Our version:\t%s s\t(%s reps).' % (t2, reps))
    return [t1,t2]

def bench_tate(E, reps=10):
    '''
    Compare our tate pairing with the tate_pairing function from sage.
    '''
    E_Fp2 = E.to_gfp2()
    EE = EllipticCurve(E_Fp2.field, [0,E_Fp2.A,0,1,0])
    p = E_Fp2.field.characteristic()
    P =  E.point_of_order(N=True, n=0, twist=False, deterministic=False)
    PP = P.get_coordinates(EE)
    Q = E_Fp2.point_of_order(N=True, n=0, twist=False, deterministic=False)
    while Q.x/Q.z in E.field:
        Q = E_Fp2.point_of_order(N=True, n=0, twist=False, deterministic=False)
    QQ = Q.get_coordinates(EE)
    t1 = cputime()
    for i in range(reps):
        PP.tate_pairing(QQ, s.N, k=2)
    t1 = cputime(t1)
    print('Sage version:\t%s s\t(%s reps).' % (t1, reps))
    t2 = cputime()
    for i in range(reps):
        pairing.tate(PP, QQ, E)
    t2 = cputime(t2)
    print('Our version:\t%s s\t(%s reps).' % (t2, reps))
    return [t1,t2]


def bench_exponentiation(setup, reps=10):
    F = setup.E0.to_gfp2().field
    t1 = cputime()
    for i in range(reps):
        (F.random_element())**((F.characteristic()**2-1)//setup.N)
    t1 = cputime(t1)
    print('Sage version:\t%s s\t(%s reps).' % (t1, reps))
    t2 = cputime()
    for i in range(reps):
        pairing.exponentiation(setup, F.random_element())
    t2 = cputime(t2)
    print('Our version:\t%s s\t(%s reps).' % (t2, reps))
    return [t1,t2]


if __name__ == '__main__':
    for key, s in setup.SETUPS.items():
        print('\n==== %s ====' % key)
        # Take a curve defined over Fp2
        P = s.E0.elligator(2)
        E1, _ = s.E0.isogeny_forward([P])
        print('E1=', E1)
        print('Benching square root function')
        bench_sqrt(s.E0)
        print('Benching tate pairing')
        bench_tate(E1)
        print('Benching pairing final exponentiation')
        bench_exponentiation(s)

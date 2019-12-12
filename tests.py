#!/usr/bin/env sage
# -*- coding: utf-8 -*-

import setup, curve, point, vdf

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

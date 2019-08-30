#-*- coding: utf-8 -*-
import time
from sage.all import *
import curve
import point

# We implement the optimal strategy when computing the isogeny walk for a toy curve

n = 12
m = 13
p = 2**n * 3**m - 1
Fp = GF(p)
Fpx = Fp['x']
x = Fpx.gen()
Fp2 = Fp.extension(x**2+1, 'u')
u = Fp2.gen()
a, b = 3151909130, 4925574689
E = EllipticCurve(Fp2, [a,b]); assert E.is_supersingular()
alpha = E.division_polynomial(2).roots()[0][0]
s = 1/sqrt(Fp(3*alpha**2 + E.a4()))
c = curve.Curve(3**m, n, 1, Integers()(3*alpha*s), 1, 0)
# change Delta for vdf use

# TEST OF CONVERSION BETWEEN WEIERSTRASS AND MONTGOMERY
def montgomery(P, alpha, s) :
    '''
    INPUT:
        * the point on the weierstrass model
        * alpha for the conversion (wiki)
        * s for the conversion (wiki)
    OUTPUT:
        * the point on the montgomery model
    '''
    x = P[0]/P[2]
    return point.Point(s*(x-alpha), 1, curve.Curve(3**m, n, 1, Integers()(3*alpha*s), 1, 0)
)

P = E.random_point()
S = E(4403911356, 2265545646) # S.order() == 4
phi = E.isogeny(S)
print 'Isogeny computation on Sage EllipticCurve obejct gives an isogeny to a curve with j-invariant j =', phi.codomain().j_invariant()

Pmg = point.Point(montgomery(P, alpha, s).x, 1, c)
Smg = montgomery(S, alpha, s)
[phiPmg] = Smg.isogeny_degree4([Pmg])
print 'Isogeny computaiton on our code gives on isogeny to a curve with j-invariant j =', phiPmg.curve.weierstrass().j_invariant()

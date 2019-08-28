#-*- coding: utf-8 -*-
import time
from sage.all import *
import curve

# We implement the optimal strategy when computing the isogeny walk
# Current status:
#  * OK for toy sage implem
#  * nOK for our code implem (not even tried)

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
A, M = c.a, s
# A, M are the coefficients of the montgomery curve

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
    return [s*(x-alpha), 1]

P = E.random_point()
P43 = (4**3) * P

import point
Pmong = point.Point(montgomery(P, alpha, s)[0], 1, c)
P43mong = (4**3)*Pmong
assert P43mong.compareXWithWeierstrass(P43)

# TEST OF VALIDITY OF isogeny_degree4 AND isogeny_degree4k
P = E(2691078013*u + 5411211369, 5875184984*u + 3991903998)
Pmong = point.Point(montgomery(P, alpha, s)[0], 1, c)

# TEST FOR isogeny_degree4
S4 = E(5212767203*u + 1541898534, 1982433626*u + 2069840204)
S4mong = point.Point(montgomery(S4, alpha, s)[0], 1, c)
[ImagePmong] = S4mong.isogeny_degree4([Pmong])
assert ImagePmong.compareXWithWeierstrass(E.isogeny(S4)(P))

# TEST FOR isogeny_degree4k_strategy
S = E(2527535319*u + 1865797195, 5803316480*u + 321180210)
Smong = point.Point(montgomery(S, alpha, s)[0], 1, c)
[phiPmong, phiP4kmong, listofcurves] = Smong.isogeny_degree4k_strategy(Pmong, 4, 'kernel4', [2,1,1])
assert phiPmong.compareXWithWeierstrass(E.isogeny(S)(P))

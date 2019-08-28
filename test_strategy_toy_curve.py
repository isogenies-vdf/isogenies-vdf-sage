#-*- coding: utf-8 -*-
import time
from sage.all import *
import curve
import point

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

Pmg = point.Point(montgomery(P, alpha, s)[0], 1, c)
P43mg = (4**3)*Pmg
assert P43mg.compareXWithWeierstrass(P43)

# TEST OF VALIDITY OF isogeny_degree4 AND isogeny_degree4k
P = E(2691078013*u + 5411211369, 5875184984*u + 3991903998)
Pmg = point.Point(montgomery(P, alpha, s)[0], 1, c)

# TEST FOR isogeny_degree4
S4 = E(5212767203*u + 1541898534, 1982433626*u + 2069840204)
S4mg = point.Point(montgomery(S4, alpha, s)[0], 1, c)
[ImagePmg] = S4mg.isogeny_degree4([Pmg])
assert ImagePmg.compareXWithWeierstrass(E.isogeny(S4)(P))

# TEST FOR isogeny_degree4k_strategy
def strategy(n, p, q):
    '''
    INPUT:
    * n the height of the tree
    * p the cost of one multiplication step
    * q the cost of one isogeny step
    OUTPUT:
    * a list corresponding to a strategy
    REMARK:
    from Luca De Feo's answer on crypto.stackexchange.com
    '''
    S = { 1: [] }
    C = { 1: 0 }
    for i in range(2, n+2):
        b, cost = min(((b, C[i-b] + C[b] + b*p + (i-b)*q) for b in range(1,i)), key=lambda t: t[1])
        S[i] = [b] + S[i-b] + S[b]
        C[i] = cost
    return S[n+1]

S = E(2527535319*u + 1865797195, 5803316480*u + 321180210)
Smg = point.Point(montgomery(S, alpha, s)[0], 1, c)
[phiPmg, phiP4kmg, listofcurves] = Smg.isogeny_degree4k_strategy(Pmg, 4, 'withoutKernel', strategy(4-1, 1, 1))
assert phiPmg.compareXWithWeierstrass(E.isogeny(S)(P))

S46 = (3**13)*E.random_point()
while S46.order() != 2**12 :
    S46 = (3**13)*E.random_point()
P = E.random_point()

S46mg = point.Point(montgomery(S46, alpha, s)[0], 1, c)
Pmg = point.Point(montgomery(P, alpha, s)[0], 1, c)

[phiPmg, phiP4kmg, listofcurves] = S46mg.isogeny_degree4k_strategy(Pmg, 6, 'kernel4', strategy(6-1, 1, 1))
print strategy(5,1,1)
assert phiPmg.compareXWithWeierstrass(E.isogeny(S46)(P))

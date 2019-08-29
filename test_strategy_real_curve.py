#-*- coding: utf-8 -*-
import time
from sage.all import *
import curve
import point

# We implement the optimal strategy when computing the isogeny walk

n = 1244
N = 0xc0256a57b1434a4970e315e3e572ad7b6b6268ca27a1bc14a5ec8d6e8f46ab63
f = 63
p = N * f * 2**n - 1
proof.arithmetic(False)
Fp = GF(p)
Fpx = Fp['x']
x = Fpx.gen()
Fp2 = Fp.extension(x**2 + 3, 'u')
u = Fp2.gen()

c = curve.Curve(63, 1244, 0xc0256a57b1434a4970e315e3e572ad7b6b6268ca27a1bc14a5ec8d6e8f46ab63, 138931309558156184106311716917677778941761847991286360325642242809534952018704195842136094062347931842162775765708572232752796610393601192925341167860358529602430304979627494497048448960083384310735203052588819895230906248500388348984991092188849520120483947949612966752973461165325952933739065855693165670941141036576698048539586409219548698834122183984266610530679658299939991747759033936995784464828547439035421618098378714023855965416127212175477937, 3, (1244-2)*28)

S = c.point_order(4**622, 2)
P = c.random_point(k=2)

t = cputime()
test1 = S.isogeny_degree4k(P, 622, 'kernel4', stop=1)
print cputime(t), ' without strat and kernel4'

t = cputime()
test2 = S.isogeny_degree4k_strategy(P, 622, 'kernel4', point.Point.strategy(622-1, 1,1), stop=1)
print cputime(t), 'with strat and kernel4'

assert test1[0] == test2[0]


t = cputime()
test1 = S.isogeny_degree4k(P, 622, 'kernel4k', stop=3)
print cputime(t), 'without strat and kernel4k'

t = cputime()
test2 = S.isogeny_degree4k_strategy(P, 622, 'kernel4k', point.Point.strategy(622-1, 1,1), stop=3)
print cputime(t), 'with strat and kernel4k'

assert test1[0] == test2[0]

#-*- coding: utf-8 -*-
import time
from sage.all import *

# We implement the optimal strategy on EllipticCurve object

from collections import deque

n = 12
m = 13
p = 2**n * 3**m - 1
Fp = GF(p)
Fpx = Fp['x']
x = Fpx.gen()
Fp2 = Fp.extension(x**2+1, 'u')
u = Fp2.gen()
E = EllipticCurve(Fp2, [3151909130, 4925574689])
assert E.is_supersingular()

# naive strategy
def naive_strat(E, S, list_of_points) :
    curve = E
    list1 = copy(list_of_points)
    R = S
    while R!=0 :
        Q = R
        while Q.order() != 4 :
            Q = 4 * Q
        phi = curve.isogeny(Q)
        R = phi(R)
        for j in range(len(list1)) :
            list1[j] = phi(list1[j])
        Q = R
        curve = phi.codomain()
    return [curve, list1]

# recursive strategy
# we do not use it, we use the iterative one
def recur_strat(E, S, list_of_points, strat) :
    if len(strat) == 0 :
        phi = E.isogeny(S)
        images = []
        for point in list_of_points :
            images.append(phi(point))
        return [phi.codomain(), images] #[phi(point) for point in list_of_points]]
    else :
        nn = strat[0]
        L = strat[1 : len(strat) + 1 - nn]
        R = strat[len(strat) + 1 - nn : len(strat)]
        T = (4**nn) * S
        [E, list1] = recur_strat(E, T, [S] + list_of_points, L)
        U = list1[0]
        [E, list2] = recur_strat(E, U, list1[1:], R)
        return [E, list2]

# optimal iterative strategy
def iter_strat(e2, E, S, list_of_points, strat) :
    QQ = deque()
    QQ.append([e2//2, S])
    i = 0
    F = E
    list1 = copy(list_of_points)
    while len(QQ) != 0 :
        [h, P] = QQ.pop()
        if h == 1 :
            phi = F.isogeny(P)
            F = phi.codomain()
            Qprime = deque()
            while len(QQ) != 0 :
                [h, P] = QQ.popleft()
                P = phi(P)
                Qprime.append([h-1, P])
            QQ = Qprime
            for j in range(len(list1)) :
                list1[j] = phi(list1[j])
        elif strat[i] > 0 and strat[i] < h :
            QQ.append([h, P])
            P = 4**(strat[i]) * P
            QQ.append([h-strat[i], P])
            i += 1
        else :
            return false
    return [F, list1]

# test
P = (3**m)*E.random_point()
while P.order() != 2**n :
    P = (3**m)*E.random_point()

S = (2**4) * P
assert S.order() == 2**8

list = [E.random_point() for i in [1,2,3]]

#optimal symmetric strategy
#  /\
# /  \
#/\  /\
strat = [2, 1, 1]
assert naive_strat(E, S, list) == recur_strat(E, S, list, strat)
assert naive_strat(E, S, list) == iter_strat(8, E, S, list, strat)

#naive peigne strategy
#  /\
# / /\
#/ / /\
strat = [3,2,1]
assert naive_strat(E, S, list) == recur_strat(E, S, list, strat)
assert naive_strat(E, S, list) == iter_strat(8, E, S, list, strat)

# We now choose fixed point and want to have the same computation with
# our code (Curve, Point, etc.)
S = E(6386625696*u + 3074559785, 3328803575*u + 3521453575)
assert S.order() == 2**n
P = E(4429946037*u + 1171750011, 4367999821*u + 5023911199)

strategy = [6,3,2,1,2,1,3,1,1,2,1]
assert naive_strat(E, S, [P]) == iter_strat(12*2, E, S, [P], strategy)
assert naive_strat(E, S, [P])[0] == E.isogeny_codomain(S)
print 'test is OK for a toy curve :-)'

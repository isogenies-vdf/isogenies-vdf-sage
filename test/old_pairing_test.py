from pairing import *

n = 64
N = 27212093
p = 2**n * N - 1
assert p.is_prime()
Fp = GF(p)
Fp2.<u> = GF(p^2, modulus = x^2 +1)

a = (189824125799020563614613076*u+164000141261416267432480249)
b = (447420608670769313888476000*u+177621199084255076097252861)
assert p == 501974515280983173562892287

assert -2*p - 1 == -1003949030561966347125784575

E = EllipticCurve(Fp2, [a,b])
assert E.is_supersingular()
P = E(267365094384829161351609777*u + 140209509831916596515193036, 487302866606018713317868073*u + 241375705642521653721213127, 1)
Q = E(348895956007647557752771548*u + 492126981510197162205833617, 254845727466868653448311619*u + 366170852225519382979594217, 1)

R = E.random_point()

l = P._line_(P, Q)
[add, [n,d]] = eval_line(P, P, Q)

assert add == 2*P and n/d == l

#tate
mi = P._miller_(Q, N)
[S, [n,d]] = miller(P, Q, N, denominator=True)

assert mi == n/d and S == N * P

"""P = (2**n) * E.random_point()
while P.order() != N :
    P = (2**n)*E.random_point()

Fp2 = GF(p**2)

EFp2 = E.change_ring(Fp2)

P = EFp2(P)

Q = (2**(2*n))*EFp2.random_point()
while Q.order() != N :
    Q = (2**(2*n))*EFp2.random_point()

R = EFp2.random_point()

l = P._line_(R, Q)
[add, [n, d]] = eval_line(R, P, Q)

assert add == P+R and n/d == l

l = P._line_(P, Q)
[add, [n,d]] = eval_line(P, P, Q)

assert add == 2*P and n/d == l


mi = P._miller_(Q, N)
[S, [n,d]] = miller(P, Q, N, denominator=True)

assert mi == n/d and S == N * P
"""

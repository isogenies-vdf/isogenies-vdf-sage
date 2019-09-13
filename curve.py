# -*- coding: utf-8 -*- 
from sage.all_cmdline import *
import point

proof.arithmetic(False)
class Curve:
    def __init__(self, f, n, N, a, ext, Delta, strategy):
        self.f = f
        self.n = n
        self.N = N
        self.p = 2**n*N*f - 1
        self.Fp = GF(self.p)
        self.Fpx = self.Fp['x']; x = self.Fpx.gen()
        self.Fp2 = self.Fp.extension(x**2+ext, 'u')
        self.a = self.Fp2(a)
        self.Delta = Delta
        self.cof_P = (self.p+1)//self.N
        self.strategy = strategy

    def __repr__(self):
        return 'Montgomery curve defined by y^2 = x^3 + (' + repr(self.a) + ')*x^2 + x over ' + repr(self.Fp2)

    __str__ = __repr__

    def random_point(self, k=1, twist=False) :
        '''
        INPUT:
        * k = 1 or 2 if you want a point defined over Fp or Fp2
        OUTPUT:
        * [x,1] a random point on the montgomery curve, in the xz-model
        '''
        F = self.Fp2 if k==2 else self.Fp
        x = F.random_element()
        #strange trick for generalize with twist
        while not(F(x**3+self.a*x**2+x).is_square()) != twist :
            x = F.random_element()
        #self is not the curve if twist=True, but anyway
        return point.Point(x, F(1), self)

    def weierstrass(self) :
        '''
        INPUT:

        OUTPUT:
        * E the elliptic curve in Weierstrass model (y^2 = x^3+a*x+b)
        '''
        return EllipticCurve(self.Fp2, [1-(self.a**2)/3, self.a*(2*(self.a**2)/9-1)/3])

    def point_order(self, k, extension_field = 2, twist = False) :
        if k == 2**(valuation(k,2)) :
            l = valuation(k,2)
            if (l > self.n) :
                print("there is no point of order 2^%d over Fp^%d" % (l, extension_field))
                return False
            cof = (self.p+1) // (2**l)
            
            if extension_field == 1 :
                # As  E(Fp) \simeq ZZ / ((p+1)/2) ZZ \times ZZ / 2 ZZ
                # there is no point of order 2^curve.n
                if l == self.n :
                    print('impossible to get a point of order 2^%d over Fp' % l)
                    return False
                else :
                    # we need to divide by 2 the cofactor because of the ZZ / 2 ZZ part above
                    cof = cof // 2
            R = cof * self.random_point(extension_field, twist)
            while not(R.is_order(2**l)) :
                R = cof * self.random_point(extension_field, twist)
        elif k == self.N :
            R = self.cof_P * self.random_point(extension_field, twist)
            while R.z == 0 :
                R = self.cof_P * self.random_point(extension_field, twist)
        return R

    def getPointFromWeierstrass(self, P) :
        '''
        INPUT:
        * P a point on the Weierstrass curve
        OUTPUT:
        * the corresponding point on the montgomery curve, in the xz-model
        '''
        if P[2] == 0 :
            return point.Point(1, 0, self)
        if P[2] == 1 :
            return point.Point(P[0] - self.a/3, 1, self)
        X, Y, Z = P[0]/P[2], P[1]/P[2], 1
        return self.getPointFromWeierstrass([X,Y,Z])

    @staticmethod
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

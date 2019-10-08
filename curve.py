# -*- coding: utf-8 -*- 
import point
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve

class Curve:
    def __init__(self, f, n, N, a, ext):
        self.f = f
        self.n = n
        self.N = N
        self.p = 2**n*N*f - 1
        self.Fp = GF(self.p)
        self.Fpx = self.Fp['x']; x = self.Fpx.gen()
        self.Fp2 = self.Fp.extension(x**2+ext, 'u')
        self.a = self.Fp2(a)
        self.cof_P = (self.p+1)//self.N

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

    def power_of_2_order_random_point(self, k, extension_field = 2, twist = False) :
        '''
            INPUT:
            * k an integer
            * extension_field an integer
            * twist a boolean
            OUTPUT:
            * R a point of the curve or its twist, of order 2**k, defined over \F_{p^k}, given in the Montgomery model
        '''
        if (k > self.n) :
            raise RuntimeError('there is no point of order 2^%d over Fp^%d' % (k, extension_field))
        cof = (self.p+1) // (2**k)

        if extension_field == 1 :
            # E(Fp) \simeq ZZ / ((p+1)/2) ZZ \times ZZ / 2 ZZ
            # there is no point of order 2^curve.n
            if k == self.n :
                raise RuntimeError('impossible to get a point of order 2^%d over Fp' % k)
            else :
                # we need to divide by 2 the cofactor because of the ZZ / 2 ZZ part above
                cof = cof // 2
        R = cof * self.random_point(extension_field, twist)
        while not(R.is_power_of_2_order_point(k)) :
            R = cof * self.random_point(extension_field, twist)
        return R

    def pairing_group_random_point(self, extension_degree = 2, twist = False) :
        '''
        INPUT:
        * extension_degree an integer
        * twist a boolean
        OUTPUT:
        * a point of order self.N defined over self.Fp^extension_degree defined over the curve or its twist.
        '''
        R = self.cof_P * self.random_point(extension_degree, twist)
        while R.z == 0 :
            R = self.cof_P * self.random_point(extension_degree, twist)
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


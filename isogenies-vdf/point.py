# -*- coding: utf-8 -*-
import sage.all
from sage.rings.integer_ring import ZZ
from collections import deque

class Point:
    def __init__(self, x, z, c):
        self.x = x
        self.z = z
        self.curve = c

    def __repr__(self):
        return '(%r:%r)' % (self.x, self.z)

    __str__ = __repr__
    

    def __eq__(self,other) :
        return self.x * other.z == self.z * other.x

    def is_zero(self):
        return self.z == 0
    
    def normalize(self) :
        '''
        INPUT:

        OUTPUT:
        * [x/z, 1] the normalized point representing P
        '''
        if self.z == 0 :
            return Point(1, 0, self.curve)
        return Point(self.x/self.z, 1, self.curve)

    def dbl(self):
        '''
        INPUT:

        OUTPUT:
        * R = [2]P, in the xz-model
        '''
        #eprint 2017/212 algo 2
        x = self.x
        z = self.z
        v1 = x+z
        v1 = v1**2
        v2 = x-z
        v2 = v2**2
        xR = v1*v2
        v1 = v1-v2
        v3 = ((self.curve.A+2)/4)*v1
        v3 = v3+v2
        zR = v1*v3
        return Point(xR, zR, self.curve)

    def add(self, Q, PmQ) :
        '''
        INPUT:
        * Q a point of the curve
        * PmQ = P-Q
        OUTPUT:
        * R = P+Q, in the xz-model
        '''
        xP, zP = self.x, self.z
        xQ, zQ = Q.x, Q.z
        xm, zm = PmQ.x, PmQ.z
        v0 = xP + zP
        v1 = xQ - zQ
        v1 = v1 * v0
        v0 = xP - zP
        v2 = xQ + zQ
        v2 = v2*v0
        v3 = v1+v2
        v3 = v3**2
        v4 = v1-v2
        v4 = v4**2
        xR = zm *v3
        zR = xm * v4
        return Point(xR, zR, self.curve)

    def __rmul__(self, k) :
        '''
        INPUT:
        * k an integer
        OUTPUT:
        * R = [k]P, in the xz-model
        '''
        if k == 0 :
            return Point(1, 0, self.curve) #Smith notation
        k = abs(k) # function does not care about sign(k)
        R0 = self
        R1 = R0.dbl()
        i=0
        R1mR0 = R0
        k_bits = ZZ(k).bits()[::-1]
        for i in range(1, len(k_bits)) :
            if k_bits[i] == 0 :
                [R0, R1] = [R0.dbl(), R0.add(R1, R1mR0)]
            else :
                [R0, R1] = [R0.add(R1, R1mR0), R1.dbl()]
        if R0.z == 0 :
            return Point(1, 0, self.curve)
        return R0

    def trace(self, image):
        '''
        Compute the trace of this point, i.e. P + π(P)

        INPUT:

        `image` is the image curve (defined over GF(p)) of the trace map

        OUTPUT:

        The trace of this point in `image`
        '''
        assert self.curve.alpha == image.alpha or self.curve.alpha == 1/image.alpha, 'Trace image curve does not match.'
        
        # using affine coordinates, maybe not the best
        x = self.x / self.z
        GFp = image.field
        # deal with the special cases self ∈ ker(π-1) ∪ ker(π+1)
        if x in GFp:
            if self in image.to_gfp2():
                # double
                tmp = Point(x, 1, image)
                return 2*tmp
            else:
                # zero point
                return Point(1, 0, image)

        else:
            # using addition formula:
            #
            #   (x,y) + (x^p,y^p) = (y^p-y)²/(x^p-x)² - A - (x^p+x)
            #
            # the code assumes GF(p²) is represented as GF(p)[i]
            y2 = (x**2 + self.curve.A*x + 1)*x

            # (y^p - y)² = (y²)^p + y² -2y^(p+1) = 2 ( Re(y²) - Norm(y) )
            normy = GFp( y2**((x.parent().characteristic() + 1) // 2) )
            num = 2 * (y2.polynomial()[0] - normy)

            # (x^p - x)² = -4 Im(x)
            den = -4 * x.polynomial()[1]**2

            return Point(num/den - image.A - 2*x.polynomial()[0], 1, image)
    
    #### Conversion to XYZ coordinates
    
    def get_coordinates(self, E):
        '''
        return the point on the EllipticCurve object (with y with a sign)
        '''
        x, z = self.x, self.z
        t = x**3/z + E.a2() * x**2 + x * z
        return E(x, self.curve._sqrt(t), z)

    #### Conversion to Weierstrass
    
    def weierstrass(self, curve, twist=False):
        '''
        INPUT:

        OUTPUT:
        * the point on the short Weierstrass curve corresponding to the montgomery curve
        '''
        if curve.is_on_gfp:
            E = curve.weierstrass(twist=twist)
        else :
            E = curve.to_gfp2().weierstrass(twist=twist)
        if self.is_zero() :
            return E(0)
        if twist :
            return E.lift_x(-self.x/self.z)
        else :
            return E.lift_x(self.x/self.z)
        # xn = self.x / self.z
        # x_w = xn + self.curve.A / 3
        # return self.curve.weierstrass().lift_x(x_w, extend=True)

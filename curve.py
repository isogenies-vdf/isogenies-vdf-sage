# -*- coding: utf-8 -*- 
import point
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from random import choice

class Curve:
    '''
    An elliptic curve in semi-Montogmery (sM) form.

    An elliptic curve in sM form is defined by the equation

        y² = x³ - (α + 1/α)x² + 1

    with α ≠ 0. It is equivalent to a Montogmery curve with coefficient
    
        A = - α - 1/α

    Consequently, α and 1/α are the absissas of two points of order 2.
    '''
    
    def __init__(self, alpha, setup):
        if (alpha + 1)**4 == 4:
            raise ValueError('Curve is singular (alpha=%d)' % alpha)
        if alpha**2 == -1:
            raise ValueError('Trying to build curve with j=1728.\n Either you did something wrong, or you are very unlucky.')
        self.alpha = alpha
        
        self.field = alpha.parent()
        self.is_on_gfp = self.field.is_prime_field()
        self.non_square = -1 if self.is_on_gfp else setup.non_square_GFp2
        self.setup = setup
        
    @property
    def A(self):
        return - self.alpha - 1/self.alpha

    @property
    def j(self):
        return 256 * (self.A**2 - 3)**3 / (self.A**2 - 4)

    def __repr__(self):
        return 'Montgomery curve y² = x³ + (%r) x² + x' % self.A

    __str__ = __repr__

    def to_gfp2(self):
        'Extend scalars to GF(p²)'
        return Curve(self.setup.GFp2(self.alpha), self.setup)
    
    def __contains__(self, P, twist=False) :
        '''
        Test whether a point `P` is on this curve.

        INPUT:
        - `P` the point to test membership for

        OPTIONS:
        - `twist` test on the twist instead

        OUTPUT:
        * True/False
        '''
        if P.z == 0:
            return True
        xz = self.field(P.x * P.z)
        xz2 = self.field(P.x + P.z)**2
        return (xz * (xz2 + (self.A - 2)*xz)).is_square() != twist

    @property
    def max_2_torsion(self):
        'The valuation of the maximal 2-Sylow'
        return self.setup.n - int(self.is_on_gfp)
    
    def elligator(self, u, twist=False):
        '''
        Elligator 2 map to the curve.

        A deterministic enconding from the base field to the curve.
        Useful for sampling points.

        The output is defined on the curve if `twist=False` (default),
        it is defined on the twist otherwise.
        '''
        if u == 0 or u == 1:
            raise ValueError('Elligator parameter must be ≥ 2.')
        ru2 = self.non_square * u**2
        ru21 = ru2 + 1
        t = self.A * ru21 * (ru2 * self.A**2 - ru21 ** 2)
        if t == 0:
            raise ValueError('Elligator paramater must not be a root of the equation')
        if t.is_square() != twist:
            return point.Point(-self.A, ru21, self)
        else:
            return point.Point(-self.A * ru2, ru21, self)

    def iterelligator(self, twist=False, deterministic=True):
        '''
        An iterator on the points of the curve, based on Elligator 2.

        If `deterministic=True` (default), it outputs Elligator(u) with 
        u enumerated according the Sage's default enumeration of the base 
        field.

        If `deterministic=False` it simply outputs random elements in the 
        image of Elligator 2 (never stops).

        The paramer `twist` (default False) is passed down to Elligator.
        '''
        it = iter(self.field)
        while True:
            u = next(it) if deterministic else self.field.random_element()
            try:
                yield self.elligator(u, twist)
            except ValueError:
                continue
        
    def point_of_order(self, N=False, n=None, twist=False, deterministic=True):
        '''
        Return a Montgomery (X:Z) point of order N^b · 2^n, where b ∈ {0, 1}.
        
        If `N=False`, then b = 0 and a point of order 2^n is drawn. Otherwise b = 1.

        If n is larger than the maximum possible order, an error is raised.
        If `n=None`, a point of maximum possible order 2^n is drawn.

        Hence, the default is to return deterministically a point on the curve of 
        order 2^n, with n as large as possible.

        Options:

        - twist: the point must be on the twist (default False). 
          This option is ignored if the curve is defined over GF(p²).
        - deterministic: sample the point deterministically (try Elligator
          function with increasing parameter). (default True)
        '''
        it = self.iterelligator(twist, deterministic)
        
        if n is None:
            n = self.max_2_torsion
        
        while True:
            P = self.setup.f * next(it)
            NP = self.setup.N * P
            NP2 = NP
            order2 = 0
            while not NP2.is_zero():
                order2 += 1
                NP2 = 2*NP2
            if order2 < n:
                continue

            if N == False:
                return 2**(order2 - n) * NP
            else:
                P = 2**(order2 - n) * P
                if (2**n * P).is_zero():
                    continue
                return P

    def _sqrt_expo(self, u):
        '''
        Compute u^((p+1)//4)
        '''
        if u == 0:
            return 0
        else:
            # exponentation to the power (p+1)/4 = 2^(n-2) * f * N
            t = u
            for i in range(self.setup.n-2) :
                t = t**2
            return t**(self.setup.f * self.setup.N)

    def _sqrt(self, u):
        '''
        Compute the square root of u if it is in Fp2.
        '''
        if u == 0:
            return 0
        if u.polynomial().degree() == 0: # u in Fp
            t = self._sqrt_expo(u)
            if t**2 == u:
                return t
            else:
                return self.to_gfp2().field.gen() * t
        else:
            if not(u.is_square()):
                raise ValueError("No square root defined over Fp2")
            # Fp2 = Fp(i) so we can compute the root in two roots over Fp :
            [a,b] = u.polynomial().list()
            n = self._sqrt_expo(a**2+b**2)
            z = (a + n)/2
            #print((a**2+b**2).is_square())
            if not(z.is_square()) :
                z = (a - n)/2
            alpha = self._sqrt_expo(z)
            beta = b/(2*alpha)
            if alpha**2 == z:
                return choice((1,-1)) * u.parent()([alpha,beta])
            else:
                return choice((1,-1)) * u.parent()([beta,-alpha])

    def isogeny_forward(self, points, principal=True):
        '''
        Compute and evaluate the degree 2 isogeny with kernel alpha.

        This isogeny is defined by the formula

            φ(x) = x (xα - 1) / (x - α)

        INPUT:

        `points` a (possibly empty) list/tuple of points onto which the isogeny
        must be evaluated.

        OPTIONS:

        `principal=True` whether or not the α' coefficient of the image 
        curve is computed by taking the principal square root of (α² - 1).

        This option is ignored if the curve is defined over GF(p²), and a random 
        square root is returned.

        OUTPUT:
    
        - The image curve
        - A (possibly empty) tuple of evaluated points.
        '''
        rt = self._sqrt(self.alpha**2 - 1)
        alpha = 2 * self.alpha * ( self.alpha + rt ) - 1
        evals = tuple(point.Point(P.x*(P.x * self.alpha - P.z), P.z*(P.x - self.alpha * P.z), self)
                     for P in points)
        return Curve(alpha, self.setup), evals

    def isogeny_backward(self, *points):
        '''
        Evaluate the dual isogeny to `isogeny_backward`.

        This isogeny is defined by the formula

            φ(x) = (x + 1)² / 4αx

        INPUT:

        `points` a tuple of points (on the image curve of `isogeny_forward`!) 
        onto which the isogeny must be evaluated.

        OUTPUT:

        The tuple of evaluated points (on this curve).
        '''
        return tuple(point.Point((P.x + P.z)**2, 4 * self.alpha * P.x * P.z, self)
                         for P in points)
    
    def weierstrass(self) :
        '''
        INPUT:

        OUTPUT:
        * E the elliptic curve in Weierstrass model (y^2 = x^3+a*x+b)
        '''
        return EllipticCurve([0,self.A,0,1,0])
#        return EllipticCurve([1-(self.A**2)/3, self.A*(2*(self.A**2)/9-1)/3])

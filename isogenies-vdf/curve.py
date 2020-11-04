# -*- coding: utf-8 -*- 
from copy import copy
import point
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from random import choice
from collections import deque

class Curve:
    '''
    An elliptic curve in semi-Montogmery (sM) form.

    An elliptic curve in sM form is defined by the equation

        y² = x³ - (α + 1/α)x² + 1

    with α ≠ 0. It is equivalent to a Montogmery curve with coefficient
    
        A = - α - 1/α

    Consequently, α and 1/α are the absissas of two points of order 2.
    '''
    
    def __init__(self, alpha, setup, force_gfp2=False):
        if (alpha + 1)**4 == 4:
            raise ValueError('Curve is singular (alpha=%d)' % alpha)
        if alpha**2 == -1:
            raise ValueError('Trying to build curve with j=1728.\n Either you did something wrong, or you are very unlucky.')
        self.alpha = alpha
        self.setup = setup
        if not(force_gfp2):
            if((alpha+1/alpha)**(setup.p-1) == 1) :
                self.field = GF(setup.p)
            else:
                self.field = alpha.parent()
            self.is_on_gfp = self.field.is_prime_field()
        else:
            self.field = self.setup.GFp2
            self.is_on_gfp = False
        self.non_square = -1 if self.is_on_gfp else setup.non_square_GFp2
        
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
        # Extend scalars to GF(p²)
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
        if self.is_on_gfp:
            rt = self._sqrt_expo(u)
            if rt**2 == u:
                return rt
            else:
                return self.setup.GFp2.gen() * rt
        else:
            if u == 0 :
                return 0
            if not(u.is_square()):
                raise ValueError("No square root defined over Fp2")
            # Fp2 = Fp(i) so we can compute the root in two roots over Fp :
            if (u in self.setup.GFp):
                rt = self._sqrt_expo(u)
                if rt**2 == u:
                    return rt
                else:
                    return self.setup.GFp2.gen() * rt
            [a,b] = u.polynomial().list()
            n = self._sqrt_expo(a**2+b**2)
            z = (a + n)/2
            if not(z.is_square()) :
                z = (a - n)/2
            alpha = self._sqrt_expo(z)
            beta = b/(2*alpha)
            if alpha**2 == z:
                r = u.parent()([alpha,beta])
            else:
                r = u.parent()([beta,-alpha])
            return choice((1,-1)) * r

        '''
        if u == 0:
            return 0
        if principal and u.polynomial().degree() == 0: # u in Fp
            t = self._sqrt_expo(u)
            if t**2 == u:
                r = t
            else:
                r = self.to_gfp2().field.gen() * t
        else:
            if not(u.is_square()):
                raise ValueError("No square root defined over Fp2")
            # Fp2 = Fp(i) so we can compute the root in two roots over Fp :
            [a,b] = u.polynomial().list()
            n = self._sqrt_expo(a**2+b**2)
            z = (a + n)/2
            if not(z.is_square()) :
                z = (a - n)/2
            alpha = self._sqrt_expo(z)
            beta = b/(2*alpha)
            if alpha**2 == z:
                r = u.parent()([alpha,beta])
            else:
                r = u.parent()([beta,-alpha])
        if principal :
            return r
        else :
            return choice((1,-1)) * r
        '''
        
    def isogeny_forward(self, points):
        '''
        Compute and evaluate the degree 2 isogeny with kernel alpha.

        This isogeny is defined by the formula

            φ(x) = x (xα - 1) / (x - α)

        INPUT:

        `points` a (possibly empty) list/tuple of points onto which the isogeny
        must be evaluated.

        OPTIONS:
        (REMOVED)
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
        new_curve = copy(self)
        new_curve.alpha = alpha
        evals = tuple(point.Point(P.x*(P.x * self.alpha - P.z), P.z*(P.x - self.alpha * P.z), new_curve)
                      for P in points)
        return new_curve, evals

    def isogeny_backward(self, points):
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

    def large_isogeny_forward(self, points, strategy, stop=0, with_dual=False):
        '''
        INPUT:
        * self the point defining the kernel of the isogeny, of degree 4**k
        * points a list of points that we want to evaluate
        * strategy a string defining the strategy to adopt: it could be hardcoded. k = len(strategy)
        * stop an integer if we want to stop before the k-th 4-isogeny.
        OUTPUT:
        * the list of the images of the points of `points`
        * the kernel of the dual isogeny
        REMARKS:
        * self needs to be such that [4**(k-1)] self  has x-coordinate != +/- 1.
        '''
        k = len(strategy)
        l = k+1
        i = 0
        E = copy(self)
        # Find `kernel` such that 2**(l-1) * `kernel` has non-zero x-coordinate
        kernel = E.point_of_order(N=False, n=l, deterministic=False)
        K2 = kernel
        for j in range(l-1):
            K2 = 2*K2
        while(K2.x == 0) :
            kernel = E.point_of_order(N=False, n=l, deterministic=False)
            K2 = kernel
            for j in range(l-1):
                K2 = 2*K2
        self.alpha = K2.x/K2.z
        # Computes a point Q that will bring the dual kernel
        if with_dual:
            # The dual isogeny has kernel phi(Q) where Q is a point of another subgroup of the `l`-torsion
            if self.is_on_gfp:
                Q = self.point_of_order(N=False, n=l, twist=True, deterministic=False)
            else :
                Q = self.point_of_order(N=False, n=l, deterministic=False)
            Q2 = Q
            for j in range(l-1):
                Q2 = 2*Q2
            while (K2.x * Q2.z == Q2.x * K2.z or K2.x * Q2.z == - Q2.x * K2.z) :
                if self.is_on_gfp:
                    Q = self.point_of_order(N=False, n=l, twist=True, deterministic=False)
                else :
                    Q = self.point_of_order(N=False, n=l, deterministic=False)
                Q2 = Q
                for j in range(l-1):
                    Q2 = 2*Q2
        images = points
        if with_dual:
            images = (Q,) + images
            list_of_curves = [self]
        queue1 = deque()
        queue1.append([l, kernel])
        while len(queue1) != 0 and l > stop :
            [h, P] = queue1.pop()
            if h == 1 :
                queue2 = deque()
                E.alpha = P.x/P.z
                while len(queue1) != 0 :
                    [h, Q] = queue1.popleft()
                    Enew, (fQ,) = E.isogeny_forward((Q,))
                    queue2.append([h-1, fQ])
                queue1 = queue2
                Enew, images = E.isogeny_forward(images)
                l -=  1
                E = Enew
                if with_dual:
                    list_of_curves.append(E)
            elif strategy[i] > 0 and strategy[i] < h :
                queue1.append([h, P])
                PP = copy(P)
                for j in range(strategy[i]):
                    PP = 2*PP
                queue1.append([h-strategy[i], PP])
                i += 1
            else :
                raise RuntimeError('There is a problem in the isogeny computation.')
        if with_dual:
            return E, images, kernel, list_of_curves[:-1]
        else :
            return E, images, kernel

    def large_isogeny_backward(self, kernel, points, strategy, listOfCurves):
        k = len(strategy)
        l = k+1
        i = 0
        images = points
        queue1 = deque()
        queue1.append([l, kernel])
        while len(queue1) != 0 :
            [h, P] = queue1.pop()
            if h == 1 :
                queue2 = deque()
                while len(queue1) != 0 :
                    [h, Q] = queue1.popleft()
                    (fQ,) = listOfCurves[l-1].isogeny_backward((Q,))
                    queue2.append([h-1, fQ])
                queue1 = queue2
                images = listOfCurves[l-1].isogeny_backward(images)
                l -=  1
            elif strategy[i] > 0 and strategy[i] < h :
                queue1.append([h, P])
                PP = copy(P)
                for j in range(strategy[i]):
                    PP = 2*PP
                queue1.append([h-strategy[i], PP])
                i += 1
            else :
                raise RuntimeError('There is a problem in the isogeny computation.')
        return images

   
    def weierstrass(self, twist=False) :
        '''
        INPUT:

        OUTPUT:
        * E the elliptic curve in Weierstrass model (y^2 = x^3+a*x+b)
        '''
        if twist and self.is_on_gfp:
            return EllipticCurve([0,-self.A,0,1,0])
        return EllipticCurve([0,self.A,0,1,0])
#        return EllipticCurve([1-(self.A**2)/3, self.A*(2*(self.A**2)/9-1)/3])

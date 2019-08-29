# -*- coding: utf-8 -*-
from sage.all import *
import curve
from collections import deque
from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism

class Point:
    def __init__(self, x, z, c):
        self.x = x
        self.z = z
        self.curve = c

    def __str__(self):
        return '[' + str(self.x) + ', ' + str(self.z) + ']'

    def __repr__(self):
        return '[' + repr(self.x) + ', ' +  repr(self.z) + ']'

    def compareXWithWeierstrass(self, other) :
        '''
        INPUT:
        * other a point of an elliptic curve of SageMath object EllipticCurve
        OUTPUT:
        * True of False if self corresponds to a point on a Weierstrass curve with the same x-coordinate as other
        '''
        field = self.curve.Fp2
        selfW = self.weierstrass()
        if selfW.curve() != other.curve() :
            # we are in an isomorphic curve
            # let's move to the right one
            a1, b1 = selfW.curve().a4(), selfW.curve().a6()
            C2 = EllipticCurve(field, [field(other.curve().a4().polynomial().list()), field(other.curve().a6().polynomial().list())])
            iso = WeierstrassIsomorphism(E=selfW.curve(), F=C2)
            selfW = iso(selfW)
        x = field(selfW[0].polynomial().list())
        z = field(selfW[2].polynomial().list())
        X = field(other[0].polynomial().list())
        Z = field(other[2].polynomial().list())
        return x * Z == z * X

    def __eq__(self,other) :
        return self.x * other.z == self.z * other.x

    def normalize(self) :
        '''
        INPUT:

        OUTPUT:
        * [x/z, 1] the normalized point representing P
        '''
        if self.z == 0 :
            return Point(1, 0, self.curve)
        return Point(self.x/self.z, 1, self.curve)

    def in_curve(self, twist=False) :
        '''
        INPUT:

        OUTPUT:
        * True/False wheter if P in E
        '''
        if self.z == 0 :
            return True
        x = self.x/self.z
        return (x**3 + self.curve.a*x**2 + x).is_square() + twist == 1

    def weierstrass(self) :
        '''
        INPUT:

        OUTPUT:
        * the point on the weierstrass curve corresponding to the montgomery curve defined with a
        '''
        x = self.x
        z = self.z
        if z == 0 :
            return self.curve.weierstrass()([0,1,0])
        xn = x/z
        # sage does not like finite fiels
        if not(xn in Integers()) :
            xn = self.curve.Fp2(xn.polynomial().list())
        if not((xn**3+self.curve.a*xn**2 + xn).is_square()) :
            print 'point on the twist'
        x_w = xn + self.curve.a/3
        return self.curve.weierstrass().lift_x(x_w)

    def equals(self, Q) :
        return self.x == Q.x and self.z == Q.z

    def dbl(self) :
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
        v3 = ((self.curve.a+2)/4)*v1
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

    def is_order(self, k) :
        if ZZ(k).is_prime() :
            return (k*self).z == 0 and self.z != 0
        if k == 1 :
            if self.z == 0 :
                return true
            else :
                return false
        if k == 2**(valuation(k, 2)) :
            l = valuation(k,2)
            return ((2**l) * self).z == 0 and (((2**l)//2) * self).z != 0
        if k == 4**(valuation(k, 4)) :
            l = valuation(k,4)
            return ((4**l) * self).z == 0 and (((4**l)//2) * self).z != 0
        return "not implemented"

    def get_P4(self, k) :
        '''
        INPUT:
        * self is a point of order 4**k of the curve
        * k an integer
        OUTPUT:
        * R = [4**(k-1)]P4powk  that is a point of order 4 of the curve
        '''
        R = self
        R_prec = R
        while R.z != 0 :
            R, R_prec = 4*R, R
        return R_prec

    def isogeny_degree4(self, Points) :
        '''
        INPUT:
        * self is a point of order 4 defining the isogeny of degree 4
        * Points a list of points for which we want the images
        OUTPUT:
        * list_images the list of the images of the points to evaluate
        REMARKS:
        * Case with  x-coordinate != ±1 from SIDH-spec.pdf
        * Case with x-coordinate = ±1 from Barreto et al. paper
        (does not work for curves j=0,1728)
        '''
        #assert isOrder4k(1, P4, a)
        list_images = []

        XP4 = self.normalize().x
        if XP4 != 1 and XP4 != -1 :
            aprime = 4*XP4**4 - 2
            curve_prime = copy(self.curve)
            curve_prime.a = aprime
            for R in Points :
                X, Z = R.x, R.z
                if Z != 0 :
                    X = X/Z
                    phiP_Xprime = -(X*XP4**2 + X -2*XP4) * X * (X*XP4 - 1)**2
                    phiP_Zprime = (X - XP4)**2*(2*X*XP4 - XP4**2-1)
                else :
                    phiP_Xprime = 1
                    phiP_Zprime = 0
                list_images.append(Point(phiP_Xprime, phiP_Zprime, curve_prime))
        else :
            aprime = 2*(self.curve.a+6) / (self.curve.a-2)
            curve_prime = copy(self.curve)
            curve_prime.a = aprime
            for R in Points :
                X, Z = R.x, R.z
                phiP_Xprime = (X+Z)**2 * (self.curve.a*X*Z + X**2 + Z**2)
                phiP_Zprime = (2-self.curve.a) * X * Z * (X-Z)**2
                list_images.append(Point(phiP_Xprime, phiP_Zprime, curve_prime))
        return list_images

    def dual_kernel_point(self, k) :
        '''
        INPUT:
        * self is a point of order 4**k defining the kernel of the isogeny
        * k an integer such that there is an isogeny of degree 4**k
        OUTPUT:
        * Q4k a point of order 4**k for which its image is the kernel of the dual isogeny
        REMARK:
        * See p.23 (Alice’s validation of Bob’s public key.) of https://eprint.iacr.org/2016/413.pdf for details
        '''
        P4 = self.get_P4(k)
        P2 = 2*P4
        Q_subgroup = False
        while not(Q_subgroup) :
            Q4k = self.curve.point_order(4**k)
            Q4 = Q4k.get_P4(k)
            Q2 = 2*Q4
            Q_subgroup = (P2.x * Q2.z != Q2.x * P2.z and P2.x * Q2.z != - Q2.x * P2.z)
            #condition equivalent to ePQ**(4**k) == 1 and ePQ**((4**k)//2) != 1
        return Q4k

    def change_iso_curve(self, a) :
        """
        INPUT:
        * a the montgomery a coefficient of the target curve
        OUTPUT:
        * the point P seen in the new curve
        """
        a1 = self.curve.Fp2(self.curve.a)
        if  (a**2 - 4) * (a1**2 - 3)**3 != (a1**2 - 4) * (a**2 - 3)**3:
            # i.e if  256*(a1**2 - 3)**3/(a1**2 - 4) != 256*(a**2 - 3)**3/(a**2 - 4)
            print 'the curves are not isomorphic !'
            return False
        E_a1 = EllipticCurve(self.curve.Fp2, [0,a1,0,1,0])
        E_a = EllipticCurve(self.curve.Fp2, [0,a,0,1,0])
        xP = self.x/self.z
        yP = sqrt(xP**3+a1*xP**2 + xP)

        P_ws = E_a1(xP, yP)
        iso = E_a1.isomorphism_to(E_a)
        curve_target = copy(self.curve)
        curve_target.a = a
        return Point(iso(P_ws)[0], 1, curve_target)

    def isogeny_degree4k_strategy(self, Q, k, method, strategy, stop=0) :
        '''
        INPUT:
        * self the point defining the kernel of the isogeny, of degree 4**k
        * Q a point that we want to evaluate
        * k such that the isogeny is of degree 4**k
        * method a string defining the method to use : withKernel4, withKernel4k or withoutKernel
        * strategy a list of integers representing the stragegy for browsing the tree
        OUTPUT:
        * phiQ the image of Q
        * phiQ4k a point generating the dual isogeny kernel   (if method = 'kernel4k')
        * listOfCurves the list of 4-isogenous curvesq        (if method = 'kernel4')
        REMARKS:
        * self needs to be such that [4**(k-1)] self  has x-coordinate != +/- 1.
        '''

        l = k
        phiP4k = self
        curve_prime = copy(self.curve)
        image_points = [Q] + [self]

        listOfCurves_a = []

        if method == 'kernel4k':
            Q4k = self.dual_kernel_point(k)
            #evaluate Q4k give the kernel of the dual isogeny
            image_points += [Q4k]
        Queue1 = deque()
        Queue1.append([k, self])

        PRINTCOUNTER = 0
        DEC = 0

        i = 0
        F = self.curve
        list1 = copy(image_points)
        while len(Queue1) != 0 and l > stop :
            [h, P] = Queue1.pop()
            if h == 1 :
                Queue2 = deque()
                while len(Queue1) != 0 :
                    [h, Q] = Queue1.popleft()
                    [Q] = P.isogeny_degree4([Q])
                    Queue2.append([h-1, Q])
                Queue1 = Queue2
                list1 = P.isogeny_degree4(list1)
                PRINTCOUNTER+=1
                if n(100*PRINTCOUNTER/k) > DEC :
                    DEC += 10
                    print floor(100*PRINTCOUNTER/k), '% of the bigstep'
                F = list1[0].curve
                if method == 'kernel4' :
                    listOfCurves_a.append(F)
                l -=  1
            elif strategy[i] > 0 and strategy[i] < h :
                Queue1.append([h, P])
                P = 4**(strategy[i]) * P
                Queue1.append([h-strategy[i], P])
                i += 1
            else :
                return false
        #output
        phiQ = list1[0]

        if method == 'kernel4k' :
            # the point defining the dual is the 3rd one of image_points evaluated by the 4-isogenies
            phiQ4k = list1[2]
        else :
            phiQ4k = ''

        if method != 'kernel4' :
            listOfCurves_a = ''
        return [phiQ, phiQ4k, listOfCurves_a]

    def isogeny_degree4k(self, Q, k, method, stop=0) :
        '''
        INPUT:
        * self the point defining the kernel of the isogeny, of degree 4**k
        * Q a point that we want to evaluate
        * k such that the isogeny is of degree 4**k
        * method a string defining the method to use : withKernel4, withKernel4k or withoutKernel
        OUTPUT:
        * phiQ the image of Q
        * phiQ4k a point generating the dual isogeny kernel   (if method = 'kernel4k')
        * listOfCurves the list of 4-isogenous curvesq        (if method = 'kernel4')
        REMARKS:
        * self needs to be such that [4**(k-1)] self  has x-coordinate != +/- 1.
        * Implementation is the naive one:
            /\
           / /\       (first possibility on DFJP 506.pdf figure 4)
          / / /\
         / / / /\
        '''
        l = k
        phiP4k = self
        curve_prime = copy(self.curve)
        image_points = [Q] + [self]
        
        listOfCurves_a = []
        
        if method == 'kernel4k':
            Q4k = self.dual_kernel_point(k)
            #evaluate Q4k give the kernel of the dual isogeny
            image_points += [Q4k]
            
        cpt = k
        while l > stop :
            P4 = image_points[1].get_P4(l)
            #verbose for large computations
            if (l*10)//k < cpt and ZZ(self.curve.p).nbits() > 100 :
                cpt = (l*10)//k
                print (10 - cpt)*10, '%  of the isogeny'
            image_points = P4.isogeny_degree4(image_points)
            if method == 'kernel4' :
                listOfCurves_a.append(image_points[0].curve)
            l = l-1
        
        #output
        phiQ = image_points[0]
        
        if method == 'kernel4k' :
            # the point defining the dual is the 3rd one of image_points evaluated by the 4-isogenies
            phiQ4k = image_points[2]
        else :
            phiQ4k = ''
        
        if method != 'kernel4' :
            listOfCurves_a = ''
        
        return [phiQ, phiQ4k, listOfCurves_a]

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
    

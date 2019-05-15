# -*- coding: utf-8 -*- 
from sage.all import *
import curve

class Point:
    def __init__(self, x, z, c):
        self.x = x
        self.z = z
        self.curve = c
        
    def __str__(self):
        return '[' + str(self.x) + ', ' + str(self.z) + ']'

    def __repr__(self):
        return '[' + repr(self.x) + ', ' +  repr(self.z) + ']'

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
        if ((xn)**3 + self.curve.a*(xn)**2 + xn).is_square() :
            return self.curve.weierstrass()(xn + self.curve.a/3, sqrt((xn)**3 + self.curve.a*(xn)**2 + xn))
        else :
            print 'point on the twist'
            return self.curve.weierstrass()(xn + self.curve.a/3, sqrt((xn)**3 + self.curve.a*(xn)**2 + xn))
    
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
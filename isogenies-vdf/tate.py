# -*- coding: utf-8 -*-
import sage.all

class Tate:
    def __init__(self,x):
        self.x=x

    def __repr__(self):
        return '(%r)' % (self.x)

    __str__ = __repr__
        
    def double_line_jac(self, S, P, a2):
        # compute the line tangente at S and evaluate it at P.
        # return the line evaluated and the point [2]S in jac coordinates.
        if P[2] == 0:
            raise ValueError("P point must be nonzero.")
        if S[2] == 0:
            return [[P.curve().base_field().one(), 1], S]
        if S[1] == 0:
            # this is a vertical line
            t1 = S[2]**2
            return [[P[0]*t1 - S[0], t1], (0,1,0)]
        # doubling formulas
        t1 = S[1]**2          # yS²
        t2 = 4*S[0]*t1        # 4 xS yS²
        t3 = 8*t1**2          # 8 yS^4
        t4 = S[2]**2          # zS²
        t5 = 3*S[0]**2 + t4 * (2*a2*S[0] + t4) #numerator after projectivation
        Z3 = 2*S[1]*S[2]
        X3 = t5**2 - 2*t2 - a2*Z3**2
        Y3 = t5*(t2 - X3) - t3
        l_den = Z3*t4       # 2 yS zS³
        l_num = l_den*P[1] - 2*t1 + t5 * (S[0] - t4*P[0])
        return [[l_num, l_den], (X3,Y3,Z3)]

    def add_line_jac(self, S, P, Q, a2):
        if Q[2] == 0:
            raise ValueError("Q must be nonzero.")
        if P[2] == 0 or S[2] == 0:
            if P[2] == 0 and S[2] == 0:
                return [[1, 1], (P[0],P[1],1)]
            if P[2] == 0:
                t1 = S[2]**2
                return [[Q[0]*t1 - S[0], t1], S]
            if S[2] == 0:
                t1 = P[2]**2
                return [[Q[0] * t1 - P[0], t1], (P[0],P[1],1)]
        if P[0] * S[2]**2 == S[0] * P[2]**2 :
            # this is a vertical line
            if P[1]*S[2]**3 == S[1]*P[2]**3 :
                # double line
                return self.double_line_jac(P, Q, a2)
            else :
                # vertical line
                return [[Q[0] - P[0], 1], (0,1,0)]
        # addition formulas
        t1 = S[2]**2      # zS²
        t2 = S[2]*t1      # zS³
        t3 = P[0]*t1       # zS² xP
        t4 = P[1]*t2       # yP zS³
        t5 = t3 - S[0]    # zS² xP - xS
        t6 = t4 - S[1]    # yP zS³ - yS
        t7 = t5**2      # (zS² xP - xS)²
        t8 = t5*t7      # (zS² xP - xS)³
        t9 = S[0] * t7    # xS*(zS² xP - xS)²
        X3 = t6**2 - (t8+2*t9+t7*t1*a2)  # (zS²yP-yS)² - (zS²xP-xS)²(zS²xP+zS²a2+xS)
        Y3 = t6*(t9 - X3) - S[1]*t8
        Z3 = S[2]*t5      # zS*(zS² xP - xS)
        # line in the general case
        #(l_den = Z3)
        l_num = Z3 * (Q[1] - P[1]) - t6 * (Q[0] - P[0])
        return [[l_num, Z3], (X3, Y3, Z3)]

    def vertical_line_jac(self, S, Q):
        '''
        S is in jacobian coordinates, Q in affine coordinates
        return the vertical line between S and -S evaluated in Q.
        '''
        if S[2] == 0 :
            return [Q.curve().base_field().one(),1]
        t1 = S[2]**2
        return [t1*Q[0] - S[0], t1]

    def miller(self, S, Q, n, a2, denominator=True) :
        # return the miller loop f_{S, n}(Q)
        if Q.is_zero():
            raise ValueError("Q must be nonzero.")
        if n.is_zero():
            raise ValueError("n must be nonzero.")
        n_is_negative = False
        if n < 0:
            n = n.abs()
            n_is_negative = True
        t_num, t_den = 1, 1
        assert S[2] == 1
        V = (S[0], S[1], 1)
        assert Q[2] == 1
        nbin = n.bits()
        i = n.nbits() - 2
        while i > -1:
            [[ell_num, ell_den], W] = self.double_line_jac(V, Q, a2)
            t_num = ell_num * t_num**2
            if denominator :
                [v_num, v_den] = self.vertical_line_jac(W, Q)
                t_num = v_den * t_num
                t_den = v_num * ell_den * t_den**2
            V = W
            if nbin[i] == 1:
                [[ell_num, ell_den], W] = self.add_line_jac(V, S, Q, a2)
                t_num = t_num*ell_num
                if denominator :
                    [v_num, v_den] = self.vertical_line_jac(W, Q)
                    t_num = v_den * t_num
                    t_den = v_num * t_den * ell_den
                V = W
            i -= 1
        if n_is_negative :
            t_num, t_den = t_den, t_num
            V[1] = -V[1]
        return [[t_num,t_den], V]

    def tate(self, P, Q, setup, denominator=True):
        a2 = P.curve().a2()
        m1, L = self.miller(P, Q, setup.N, a2, denominator)
        m2 = m1[0]/m1[1]
        # efficient final exponentiation
        x1 = (m2**setup.f)
        for i in range(setup.n+1) :
            x1 = x1**2
        x2 = x1
        x3 = x2**(setup.f*setup.N)
        for i in range(setup.n-1) :
            x3 = x3**2
        return x3/x2

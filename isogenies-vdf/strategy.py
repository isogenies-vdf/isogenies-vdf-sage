# -*- coding: utf-8 -*-

class Strategy:
    def __init__(self, x):
        self.x = x
        
    def generate(self, n, p, q):
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
        for i in range(2, n+1):
            b, cost = min(((b, C[i-b] + C[b] + b*p + (i-b)*q) for b in range(1,i)), key=lambda t: t[1])
            S[i] = [b] + S[i-b] + S[b]
            C[i] = cost
        return S[n]


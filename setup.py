# -*- coding: utf-8 -*-

import curve
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.finite_field_constructor import FiniteField

class TrustedSetup:
    '''
    A class representing the parameters provided by a trusted setup, namely:

    - A prime p of the form (f · 2^n · N - 1)
    - A starting curve E0 of equation y² = x³ + Ax² + x
    - alpha?
    '''
    
    def __init__(self, f, n, N, alpha):
        self.f = ZZ(f)
        self.n = ZZ(n)
        self.N = ZZ(N)
        self.alpha = self.GFp(alpha)

    def __repr__(self):
        return 'Starting %s over finite field of size (2^%d·N·%d - 1), with N a %d-bits integer.' % (self.E0, self.n, self.f, self.N.nbits())

    @property
    def p(self):
        'The prime p = (f · 2^n · N - 1) defining the base field'
        return ZZ(2)**self.n * self.N * self.f - 1

    @property
    def cofactor(self):
        'The cofactor of N in the curve order, namely (f · 2^n)'
        return ZZ(2)**self.n * self.f
    
    @property
    def GFp(self):
        'The base field GF(p)'
        return FiniteField(self.p)

    @property
    def GFp2(self):
        'The degree two extension of the base field GF(p²) ≃ GF(p)[i] / (i² + 1)'
        R = self.GFp['x']
        poly = R.gen()**2 + 1
        return self.GFp.extension(poly, name='i')

    @property
    def non_square_GFp2(self):
        'A non-square element in GFp2'
        i = self.GFp2.gen()
        return next(i + x for x in self.GFp if not (i + x).is_square())
    
    @property
    def E0(self):
        'The staring curve'
        return curve.Curve(self.alpha, self)

# Some predefined setups
# The curve is always y² = x³ - 6x² + x (not safe), and α = 3 + 2√2
SETUPS = {
    'p14-toy' : TrustedSetup(f=1, n=8, N=53, alpha=5052),
    'p89-toy' : TrustedSetup(f=1, n=64, N=27212093, alpha=213532735709189898821123859),
    'p1506' : TrustedSetup(
        f = 63,
        n = 1244,
        N = 0xc0256a57b1434a4970e315e3e572ad7b6b6268ca27a1bc14a5ec8d6e8f46ab63,
        alpha = 340938674324042923738404365843495497281528715196988250806912265839837005979193060681496022381524514842692063173327721071390364939958198680950956861494810503597131815339059968476668390170096408017588488924813194150837074922948389535035013208442451816943521842939304662425358712141280829220081897275422358547644761975364052530387051257488300113871975266049432811198494556196246120225872599544980805045972355007773375688869019282622665425400703993651940608,
    ),
}
def setup_list():
    'Print the list of predefined setups'
    return '\n'.join('%s:\n  %s' % i for i in SETUPS.items())


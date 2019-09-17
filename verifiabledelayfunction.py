# -*- coding: utf-8 -*- 

class VerifiableDelayFunction:
    def __init__(self, method, strategy, curve, delay):
        self.method = method
        self.strategy = strategy
        self.curve = curve
        self.delay = delay

    def __repr__(self):
        return 'Verifiable delay function with delay ' + repr(self.delay)

    __str__ = __repr__

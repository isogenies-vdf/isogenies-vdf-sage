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

    def walk_back(self, Q, curvesPath, kernelsOfBigSteps) :
        k = len(kernelsOfBigSteps)
        T = Q
        c_t = copy(Q.curve)

        # At the moment, isogeny codomain is "up to isomorphism".
        # From the kernels given in setup, I need to move them into the 
        # right curve.
        # That is why `kernel4` is not efficient (change_iso_curve times 
        # the number of steps in the walk on the graph. In kernel4k, it 
        # is reduced by a factor 1000.
        if self.method == 'kernel4k' :
            for R in kernelsOfBigSteps :
                R = R.change_iso_curve(T.curve.a)
                [T, kernelPoint, listOfCurves] = R.isogeny_degree4k(T, k, method='withoutKernel', strategy=self.strategy)
        elif self.method == 'kernel4' :
            for c1 in curvesPath:
                R = Point(1, 1, c1).change_iso_curve(T.curve.a)
                [T] = R.isogeny_degree4([T])
        return T.change_iso_curve(self.curve.a)


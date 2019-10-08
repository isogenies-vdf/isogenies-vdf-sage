# -*- coding: utf-8 -*- 
from copy import copy

class VerifiableDelayFunction:
    def __init__(self, method, strategy, curve, delay):
        self.method = method
        self.strategy = strategy
        self.curve = curve
        self.delay = delay

    def __repr__(self):
        return 'Verifiable delay function with delay ' + repr(self.delay)

    __str__ = __repr__

    def setup_walk(self, extensionDegree, P, stop, conditionP4=True) :
        '''
        INPUT:
        * extensionDegree for the extension of the field where we choose points for the walk
        * P a point of self.curve for which we want the image
        * conditionP4 a boolean if we want to look at P4.x = +/- 1
        OUTPUT:
        * phiP the image of P by the random isogeny
        * dualKernels the list of the dual isogenies kernels
        REMARK:
        * Do not walk to the j=0,1728 curve at first step! The formulas of [4]-isogenies for these curves do not work for the moment.
        '''
        # Using isogeny_degree4k, we compute k steps in a row
        k = len(self.strategy)
        # We need to do delay // k steps
        assert self.delay % k == 0
        nbSteps = self.delay // k

        c = copy(self.curve)
        images = [P]
        dualKernels = []

        first = True

        for i in range(nbSteps) :
            if first :
                # we check the first step does not go to j=0 or 1728 curve
                j = 0
                while j == 0 or j == 1728 :
                    kernel = c.power_of_2_order_random_point(2*k, extensionDegree, False)
                    P4 = kernel.get_P4(k)
                    xP4 = P4.normalize().x
                    while (xP4 == 1 or xP4 == -1) and conditionP4 :
                        kernel = c.power_of_2_order_random_point(2*k, extensionDegree, False)
                        P4 = kernel.get_P4(k)
                        xP4 = P4.normalize().x
                    curve_onestep = P4.isogeny_degree4(images)[0].curve
                    j = curve_onestep.weierstrass().j_invariant()
                first = False
            else :
                # kernel defines an isoeny of degree 4**k
                kernel = c.power_of_2_order_random_point(2*k, extensionDegree, False)
                P4 = kernel.get_P4(k)
                xP4 = P4.normalize().x
                while (xP4 == 1 or xP4 == -1) and conditionP4 :
                    kernel = c.power_of_2_order_random_point(2*k, extensionDegree, False)
                    P4 = kernel.get_P4(k)
                    xP4 = P4.normalize().x
            # the image of dual_kernel will define the dual isogeny
            dual_kernel = kernel.dual_kernel_point(k)
            images.append(dual_kernel)
            images = kernel.isogeny_degree4k(images, self.strategy, stop)
            c = images[0].curve
            dualKernels.append(images.pop())
            phiP = images[0]
        dualKernels = dualKernels[::-1]
        return [phiP, dualKernels]

    def evaluation_walk(self, Q, dualKernels, stop):
        '''
        INPUT:
        * Q the second point of the protocol
        * dualKernels from the setup
        OUTPUT:
        * hat_phiQ the image of Q by the dual walk
        '''
        k = len(dualKernels)
        T = Q
        c_t = copy(Q.curve)

        # TODO this is not efficient.
        for R in dualKernels :
            R = R.change_iso_curve(T.curve.a)
            [T] = R.isogeny_degree4k([T], self.strategy, stop)
        return T.change_iso_curve(self.curve.a)


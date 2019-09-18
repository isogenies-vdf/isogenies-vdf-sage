# -*- coding: utf-8 -*-
import point
import curve
import argparse
from sage.rings.integer_ring import ZZ
from sage.misc.misc import cputime

parser = argparse.ArgumentParser()
parser.add_argument("--protocol", type=str, default="fp",
                    help="choose the VDF protocol to use")
parser.add_argument("--method", type=str, default="kernel4k",
                    help="choose the method to store the setup walk")
parser.add_argument("--pSize", type=str, default="p14-toy",
                    help="determine the size of the prime p to use")
parser.add_argument("--nbIterations", type=int, default=12,
                    help="set the number of iterations of the VDF")
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

args = parser.parse_args()

protocol = args.protocol
method = args.method
pSize = args.pSize
nbIterations = int(args.nbIterations)

verbose =  args.verbose

if verbose:
    print("verbosity turned on")

if pSize == 'p14-toy' :
    f = 1
    n = 8
    N = 53
    a = 10088
    alpha = 1
    strategy = [2, 1, 1]
if pSize == 'p89-toy' :
    f = 1
    n = 64
    N = 27212093
    a = 6
    alpha = 1
    strategy = [16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]
if pSize == 'p1506' :
    f = 63
    n = 1244
    N = 0xc0256a57b1434a4970e315e3e572ad7b6b6268ca27a1bc14a5ec8d6e8f46ab63
    a = 138931309558156184106311716917677778941761847991286360325642242809534952018704195842136094062347931842162775765708572232752796610393601192925341167860358529602430304979627494497048448960083384310735203052588819895230906248500388348984991092188849520120483947949612966752973461165325952933739065855693165670941141036576698048539586409219548698834122183984266610530679658299939991747759033936995784464828547439035421618098378714023855965416127212175477937
    alpha = 3
    strategy = [256, 128, 110, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 46, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 14, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 6, 4, 2, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 128, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 64, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1]

# Choice of Delta depending on nbIterations and n
Delta = nbIterations
if protocol == 'fp' :
    while Delta % (n-2) != 0 or (Delta // (n-2)) % 2 != 0:
        Delta += 1
else :
    while Delta % n != 0 or (Delta // n) % 2 != 0:
        Delta += 1
print('Delta = %d' % Delta)

c = curve.Curve(f, n, N, a, alpha)

print(c)

from fpverifiabledelayfunction import FpVerifiableDelayFunction

VDF = FpVerifiableDelayFunction(method, c, strategy, Delta)

print(VDF)


time = cputime()
setup = setup(c, verbose, method)
time = cputime(time)
print('setup timing: %.5f seconds.' % time)
[P, c_prime, curvesPath, kernelsOfBigSteps, phiP] = setup

#Generating a point Q
#NOT IN THE SAME SUBGROUP AS phiP !!!
phiP_ws = phiP.weierstrass()
e_phiP_Q = 1
while e_phiP_Q == 1 :
    Q = c_prime.cof_P * c_prime.random_point(1 if protocol=='fp' else 2)
    while Q.z == 0 :
        Q = c_prime.cof_P * c_prime.random_point(1 if protocol=='fp' else 2)
    Q_ws = Q.weierstrass()
    e_phiP_Q = Q_ws.weil_pairing(phiP_ws, ZZ(c_prime.N))

#EVAL
time = cputime()
Tr_hat_phiQ = evaluate(c, setup, Q, verbose, method)
time = cputime(time)
print('eval timing: %.5f seconds.' % time)

#VERIFY
time = cputime()
ver = verify(c, setup, Q, Tr_hat_phiQ)
time = cputime(time)
print('verif timing: %.5f seconds.' % time)
print('###############')
if ver :
    print('#verif OK  :-)#')
else :
    print('#verif nOK :-(#')
print('###############')


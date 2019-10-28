# -*- coding: utf-8 -*-
import logging
import point
import curve
import argparse
from sage.rings.integer_ring import ZZ
from sage.misc.misc import cputime
from verifiabledelayfunction import VerifiableDelayFunction
from fpverifiabledelayfunction import FpVerifiableDelayFunction
from fp2verifiabledelayfunction import Fp2VerifiableDelayFunction

parser = argparse.ArgumentParser()
parser.add_argument("--protocol", type=str, default="fp",
                    help="choose the VDF protocol to use")
parser.add_argument("--method", type=str, default="kernel4k",
                    help="choose the method to store the setup walk")
parser.add_argument("--pSize", type=str, default="p14-toy",
                    help="determine the size of the prime p to use")
parser.add_argument("--nbIterations", type=int, default=12,
                    help="set the number of iterations of the VDF")
parser.add_argument("--loglevel", type=str, default="INFO",
                    help="determines the severity for the logs")

args = parser.parse_args()

numeric_level = getattr(logging, (args.loglevel).upper(), None)
if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: %s' % args.loglevel)
#logging.basicConfig(filename='vdf.log', filemode='w', level=numeric_level, format='%(asctime)s %(message)s', datefmt='%d/%m/%Y %I:%M:%S %p')
logging.basicConfig(filename='vdf.log', filemode='w', level=numeric_level, format='%(message)s')

protocol = args.protocol
method = args.method
pSize = args.pSize
nbIterations = int(args.nbIterations)

if pSize == 'p14-toy' :
    f = 1
    n = 8
    N = 53
    a = 10088
    alpha = 1
if pSize == 'p89-toy' :
    f = 1
    n = 64
    N = 27212093
    a = 6
    alpha = 1
if pSize == 'p1506' :
    f = 63
    n = 1244
    N = 0xc0256a57b1434a4970e315e3e572ad7b6b6268ca27a1bc14a5ec8d6e8f46ab63
    a = 138931309558156184106311716917677778941761847991286360325642242809534952018704195842136094062347931842162775765708572232752796610393601192925341167860358529602430304979627494497048448960083384310735203052588819895230906248500388348984991092188849520120483947949612966752973461165325952933739065855693165670941141036576698048539586409219548698834122183984266610530679658299939991747759033936995784464828547439035421618098378714023855965416127212175477937
    alpha = 3

c = curve.Curve(f, n, N, a, alpha)

# The delay will depend on nbIterations *and n*
Delta = nbIterations

logging.info('Protocol: %s', protocol)
logging.info('Method: %s', method)
logging.info('Prime size: %s', pSize)
logging.info('Number of steps: %s', str(Delta))

if protocol == 'fp' :
    VDF = FpVerifiableDelayFunction(method, c, Delta)
else :
    VDF = Fp2VerifiableDelayFunction(method, c, Delta)

time = cputime()
SETUP = VDF.setup()
time = cputime(time)
print('setup timing: %.5f seconds.' % time)
[P, phiP, dualKernels] = SETUP

c2 = dualKernels[0].curve

logging.info('Setup:\t\t\t\t%s seconds', str(time))

#Generating a point Q
#NOT IN THE SAME SUBGROUP AS phiP !!!
phiP_ws = phiP.weierstrass()

e_phiP_Q = 1
while e_phiP_Q == 1 :
    Q = c2.cof_P * c2.random_point(1 if protocol=='fp' else 2, False)
    while Q.z == 0 :
        Q = c2.cof_P * c2.random_point(1 if protocol=='fp' else 2, False)
    Q_ws = Q.weierstrass()
    e_phiP_Q = Q_ws.weil_pairing(phiP_ws, ZZ(c2.N))

#EVAL
time = cputime()
Tr_hat_phiQ = VDF.evaluate(Q, dualKernels)
time = cputime(time)
print('eval timing: %.5f seconds.' % time)
logging.info('Evaluation:\t\t\t%s seconds', str(time))

#VERIFY
time = cputime()
ver = VDF.verify(P, phiP, Q, Tr_hat_phiQ)
time = cputime(time)
print('verif timing: %.5f seconds.' % time)
logging.info('Verification:\t\t\t%s seconds', str(time))

print('###############')
if ver :
    print('#verif OK  :-)#')
    logging.info('\t\t\t\tVerification OK')
else :
    print('#verif nOK :-(#')
    logging.info('\t\t\t\tVerification NOT OK')
    #logging.info(setup = %s', str(SETUP))
    logging.info('Tr_hat_phiQ = %s', str(Tr_hat_phiQ))

print('###############')

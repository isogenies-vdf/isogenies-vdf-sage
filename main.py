#!/usr/bin/env sage
# -*- coding: utf-8 -*-

import logging
import point
import curve
from setup import SETUPS, setup_list
import argparse, sys
from sage.rings.integer_ring import ZZ
from sage.misc.misc import cputime
from verifiabledelayfunction import VerifiableDelayFunction
from fpverifiabledelayfunction import FpVerifiableDelayFunction
from fp2verifiabledelayfunction import Fp2VerifiableDelayFunction

parser = argparse.ArgumentParser()
parser.add_argument("--setup", type=str, default="p14-toy",
                    help="choose which trusted setup to use")
parser.add_argument("-l", "--list-setups", action='store_true',
                    help="list available trusted setups")
parser.add_argument("--fp", action='store_const', const="fp", dest='protocol',
                    help="run the GF(p) variant of the vdf (default)")
parser.add_argument("--fp2", action='store_const', const="fp2", dest='protocol',
                    help="run the GF(p²) variant of the vdf")
parser.add_argument("--method", type=str, default="kernel4k",
                    help="choose the method to store the setup walk")
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

protocol = args.protocol or 'fp'
method = args.method
try:
    setup = SETUPS[args.setup]
except KeyError:
    raise ValueError('Unknown trusted setup %s. Available setups:\n',
                         args.setup, setup_list())
Delta = int(args.nbIterations)

if args.list_setups:
    print('Available setups:')
    print(setup_list())
    sys.exit(0)
else:
    logging.info('Setup: %s', setup)
    logging.info('Protocol: %s', protocol)
    logging.info('Method: %s', method)
    logging.info('Number of steps: %s', str(Delta))

    if protocol == 'fp' :
        VDF = FpVerifiableDelayFunction(setup, method, Delta)
    else :
        VDF = Fp2VerifiableDelayFunction(setup, method, Delta)

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

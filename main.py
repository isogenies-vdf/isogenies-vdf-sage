#!/usr/bin/env sage
# -*- coding: utf-8 -*-

from sage.structure.proof.all import arithmetic
arithmetic(False)

import vdf
from setup import SETUPS, setup_list
import argparse, sys, logging
from sage.misc.misc import cputime

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--setup", type=str, default="p14-toy",
                    help="choose which trusted setup to use")
parser.add_argument("-l", "--list-setups", action='store_true',
                    help="list available trusted setups")
parser.add_argument("--fp", action='store_const', const="fp", dest='protocol',
                    help="run the GF(p) variant of the vdf (default)")
parser.add_argument("--fp2", action='store_const', const="fp2", dest='protocol',
                    help="run the GF(p²) variant of the vdf")
parser.add_argument("--method", type=str, default="kernel4k",
                    help="choose the method to store the setup walk")
parser.add_argument("-T", "--nbIterations", type=int, default=12,
                    help="set the number of iterations of the VDF")
parser.add_argument("--loglevel", type=str, default="INFO",
                    help="determines the severity for the logs")

args = parser.parse_args()

numeric_level = getattr(logging, (args.loglevel).upper(), None)
if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: %s' % args.loglevel)
logging.basicConfig(level=numeric_level, format='%(message)s')

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
    logging.info('N = %r', setup.N)
    logging.info('Protocol: %s', protocol)
    logging.info('Method: %s', method)
    logging.info('Number of steps: %s\n', str(Delta))

    vdfclass = vdf.VDF_GFp if protocol == 'fp' else vdf.VDF_GFp2

    logging.info('Setting up VDF')
    time = cputime()
    vdf = vdfclass(setup, Delta)
    time = cputime(time)
    logging.info('Total time spent: %.5f s, throughput %.5f isog/ms' % (time, Delta / time / 1000))
    
    Q = vdf.random_input()

    logging.info('Evaluating')
    time = cputime()
    fQ = vdf.evaluate(Q)
    time = cputime(time)
    logging.info('Total time spent: %.5f s, throughput %.5f isog/ms' % (time, Delta / time / 1000))

    logging.info('Verifying')
    time = cputime()
    res = vdf.verify(Q, fQ)
    time = cputime(time)
    logging.info('Total time spent: %.5f s' % time)

    assert res, 'Verification failed!'

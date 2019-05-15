# -*- coding: utf-8 -*-
import time
import point
import curve

from sage.all import *

paramToSet = ['fp', 'kernel4k', 'p14-toy']

for i in range(1, len(sys.argv)):
    paramToSet[i-1] = sys.argv[i]
if paramToSet[0] != 'fp' and paramToSet[0] != 'fp2' :
    paramToSet[0] = 'fp'
if paramToSet[1] != 'kernel4' and paramToSet[1] != 'kernel4k' :
    paramToSet[1] = 'kernel4k'
if paramToSet[2] != 'p14-toy' and paramToSet[2] != 'p89-toy' and paramToSet[2] != 'p1506' :
    paramToSet[2] = 'p14-toy'

[protocol, method, pSize] = paramToSet

verbose = False

file = open("timing_" + protocol + "_" + method + "_" + pSize + ".txt", "a")

if pSize == 'p14-toy' :
    c = curve.Curve(1, 8, 53, 10088, 1, (8-2)*12)
if pSize == 'p89-toy' :
    c = curve.Curve(1, 64, 27212093, 6, 1, (64-2)*12)
if pSize == 'p1506' :
    c = curve.Curve(63, 1244, 0xc0256a57b1434a4970e315e3e572ad7b6b6268ca27a1bc14a5ec8d6e8f46ab63, 138931309558156184106311716917677778941761847991286360325642242809534952018704195842136094062347931842162775765708572232752796610393601192925341167860358529602430304979627494497048448960083384310735203052588819895230906248500388348984991092188849520120483947949612966752973461165325952933739065855693165670941141036576698048539586409219548698834122183984266610530679658299939991747759033936995784464828547439035421618098378714023855965416127212175477937, 3, (1244-2)*28)
if protocol == 'fp' :
    # for the Fp protocol, choose Delta = (n-2) * evenNumber
    assert c.Delta % (c.n - 2) == 0
    assert (c.Delta // (c.n - 2)) % 2 == 0
else :
    # for the Fp2 protocol, choose Delta = n * number (?TOCHECK)
    c.Delta = c.Delta + 2 * (c.Delta // (c.n - 2))
    assert c.Delta % c.n == 0

file.write("VDF over " + protocol + " with log_2(p) = " + str(ZZ(c.p).nbits()) + " and log_2(T) = " + str(ZZ(c.Delta).nbits()) + ".\n")
if method == 'kernel4' :
    file.write('Points of [4]-torsion stored.\n\n')
elif method == 'kernel4k' :
    file.write('Points of [4^k]-torsion stored.\n\n')

if protocol == 'fp' :
    from vdf_fp_setup import *
else :
    from vdf_fp2_setup import *

time = cputime()
setup = vdf_setup(c, verbose, method)
time = cputime(time)
print 'setup timing: ', time, ' seconds.'
[P, c_prime, curvesPath, kernelsOfBigSteps, phiP] = setup
file.write("Setup:\t" + str(time) + " seconds.\n")

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
if protocol == 'fp' :
    from vdf_fp_eval import *
else :
    from vdf_fp2_eval import *

time = cputime()
Tr_hat_phiQ = vdf_eval(c, setup, Q, verbose, method)
time = cputime(time)
print 'eval timing: ', time, ' seconds.'
file.write('Eval:\t'+ str(time) + ' seconds.\n')


#VERIFY
if protocol == 'fp' :
    from vdf_fp_verif import *
else :
    from vdf_fp2_verif import *

time = cputime()
ver = vdf_verif(c, setup, Q, Tr_hat_phiQ)
time = cputime(time)
print 'verif timing: ', time, ' seconds.'
file.write('Verif:\t' + str(time) + ' seconds.\n')
print '###############'
if ver :
    print '#verif OK  :-)#'
    file.write('verif ok\n')
else :
    print '#verif nOK :-(#'
    file.write('verif nok\n')
    file.write("setup = " + str(setup)+ "\n")
    file.write("Tr_hat_phiQ = " + str(Tr_hat_phiQ) + "\n")
print '###############'

file.write("\n")
file.close()

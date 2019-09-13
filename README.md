# Isogenies VDF Paper

https://eprint.iacr.org/2019/166.pdf

# How to compute our VDF ?
`sage vdf.py`
optional arguments:
  -h, --help            show this help message and exit
  --protocol PROTOCOL   choose the VDF protocol to use
  --method METHOD       choose the method to store the setup walk
  --pSize PSIZE         determine the size of the prime p to use
  --nbIterations NBITERATIONS
                        set the number of iterations of the VDF
  -v, --verbose         increase output verbosity

# Issues.

- How to add two points in x-model montgomery ? (in order to compute the Trace at
the end of vdf_eval step.
I fix it going to Weierstrass model and computing phi(Q) + frob(phi(Q)) and
phi(Q) - frob(phi(Q)) because maybe there is a problem of sign on y. I find 
the right one checking if the coordinates lie in Fp.
- One pairing on the verify step can be done before the answer of vdf_eval.
- Do I really need all the a-coefficients of the steps in the output of setup ?
When I compute the evaluation step, I need to switch to the right isomorphic 
curve and so need the a coefficient.
Maybe it can be recomputed during eval... (?)

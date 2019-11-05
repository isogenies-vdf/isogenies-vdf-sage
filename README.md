# Isogenies VDF Paper

https://eprint.iacr.org/2019/166.pdf

# How to compute our VDF ?
`sage -python3 vdf.py --protocol PROTOCOL --method METHOD --pSize PSIZE --nbIterations NBITERATIONS --loglevel LOGLEVEL`
where:
- `PROTOCOL` is `fp` or `fp2`
- `METHOD` is `kernel4k` or `kernel4`
- `PSIZE` is `p14-toy`, `p89-toy` or `p1506`
- `NBITERATIONS` is the number of 2-isogenies to compute for the evaluation.
- `LOGLEVEL` is the severity for the logging.

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

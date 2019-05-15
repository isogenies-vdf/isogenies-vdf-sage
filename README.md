# How to compute our VDF ?
`sage vdf.py protocol method pSize` where:
- `protocol` is `fp` or `fp2`,
- `method` is `kernel4` or `kernel4k`
- `pSize` is `p14-toy`, `p89-toy` or `p1506`.

# Issues.

- How to add two points in x-model montgomery ? (in order to compute the Trace at
the end of vdf_eval step.
I fix it going to Weierstrass model and computing phi(Q) + frob(phi(Q)) and
phi(Q) - frob(phi(Q)) because maybe there is a problem of sign on y. I find 
the right one checking if the coordinates lie in Fp.
- Ate pairing is not efficent for one of the pairings in the Fp2 case (we need 
to compute the denominators...).
- Optimal strategy for computing a [4^k]-isogeny is not done.
- One pairing on the verify step can be done before the answer of vdf_eval.
- Do I really need all the a-coefficients of the steps in the output of setup ?
When I compute the evaluation step, I need to switch to the right isomorphic 
curve and so need the a coefficient.
Maybe it can be recomputed during eval... (?)

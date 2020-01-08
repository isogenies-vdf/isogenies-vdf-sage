# Isogenies VDF Paper

https://eprint.iacr.org/2019/166.pdf

# How to compute our VDF ?
`sage main.py -s PSIZE --PROTOCOL -T NBITERATIONS`
where:
- `PROTOCOL` is `fp` or `fp2`
- `PSIZE` is `p14-toy`, `p89-toy` or `p1506`
- `NBITERATIONS` is the number of 2-isogenies to compute for the evaluation.

# Issues.

- One pairing on the verify step can be done before the answer of vdf_eval.
- Do I really need all the a-coefficients of the steps in the output of setup ?
When I compute the evaluation step, I need to switch to the right isomorphic 
curve and so need the a coefficient.
Maybe it can be recomputed during eval... (?)

# Isogenies VDF Paper

https://eprint.iacr.org/2019/166.pdf

# How to compute our VDF ?
`sage main.py -s PSIZE --PROTOCOL -T NBITERATIONS`
where:
- `PROTOCOL` is `fp` or `fp2`
- `PSIZE` is `p14-toy`, `p89-toy` or `p1506`
- `NBITERATIONS` is the number of 2-isogenies to compute for the evaluation.

# Test
`sage test/test_curve.py`
`sage test/test_point.py`
`sage test/test_tate.py`
`sage test/test_verifiable_delay_function.py`

# Bench
`sage bench/benchs.py` (not working for the moment)

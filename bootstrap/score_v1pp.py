#!/usr/bin/env python3
"""Score V1'' exactly as registered in SPEC.md. Run from repo root:
    python bootstrap/score_v1pp.py > bootstrap/V1pp.txt"""
import itertools
import json
import os

import numpy as np
from scipy.stats import chi2, norm

HERE = os.path.dirname(os.path.abspath(__file__))
SEEDS = (101, 202, 303)
ok_all = True
for dim in ('2d', '4d'):
    R = [json.load(open(os.path.join(HERE, 'v1pp', f'{dim}_{s}.json'))) for s in SEEDS]
    N = R[0]['N']
    Z = np.array([r['z'] for r in R])
    zc = Z.sum(0) / np.sqrt(len(SEEDS))
    zN = norm.ppf(1 - 0.005 / N)
    a = np.max(np.abs(zc)) < zN
    print(f'{dim.upper()}: {N} equations, 3 seeds')
    print(f'  (a) combined max |z| = {np.max(np.abs(zc)):.2f}  vs  z_N = {zN:.2f}  -> {a}'
          f'   (worst eq {int(np.argmax(np.abs(zc)))}: per-seed z = '
          + ', '.join(f'{x:+.2f}' for x in Z[:, int(np.argmax(np.abs(zc)))]) + ')')
    ok = a
    if dim == '2d':
        lim = 3 / np.sqrt(N)
        for i, j in itertools.combinations(range(3), 2):
            r = np.corrcoef(Z[i], Z[j])[0, 1]
            b = abs(r) < lim
            ok &= b
            print(f'  (b) corr(seed {SEEDS[i]}, seed {SEEDS[j]}) = {r:+.3f}  (limit +-{lim:.3f}) -> {b}')
    else:
        res = np.mean([r['res'] for r in R], axis=0)
        covs = []
        for r in R:
            reps = np.array(r['reps'])
            nb = len(reps)
            dlt = reps - reps.mean(0)
            covs.append((nb - 1) / nb * dlt.T @ dlt)
        cov = np.mean(covs, axis=0) / len(R)
        c2 = float(res @ np.linalg.solve(cov, res))
        p = 1 - chi2.cdf(c2, N)
        c = p > 0.001
        ok &= c
        print(f'  (c) full-covariance chi2 = {c2:.1f} / {N}, p = {p:.3f} -> {c}')
    d = all(abs(r['control_z']) > 5 for r in R)
    ok &= d
    print(f'  (d) negative control z = ' + ', '.join(f"{r['control_z']:+.1f}" for r in R) + f' -> {d}')
    print(f'  V1\'\' {dim.upper()} -> {"PASS" if ok else "FAIL"}\n')
    ok_all &= ok
print(f"V1'' OVERALL -> {'PASS' if ok_all else 'FAIL: equations treated as wrong; no bound reported'}")

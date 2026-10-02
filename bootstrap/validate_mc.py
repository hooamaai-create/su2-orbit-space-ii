#!/usr/bin/env python3
"""BOOT-01 gate V1: every generated loop equation must hold on Monte Carlo data.

Run from repo root:  python bootstrap/validate_mc.py > bootstrap/V1.txt
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'contrast'))
from loops import Algebra, nonbacktracking_paths, endpoint  # noqa: E402
from su2_lattice import Lattice, qdag, qmul  # noqa: E402

D, L, BETA, THERM, NCFG, SKIP, NB = 4, 6, 2.3, 200, 120, 5, 12
LENGTHS = (4, 6)
SEED = int(os.environ.get('BOOT_SEED', 11))
if len(sys.argv) > 1:                       # diagnostic overrides: D L beta ncfg lengths...
    D, L, BETA, NCFG = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    LENGTHS = tuple(int(x) for x in sys.argv[5:]) or LENGTHS
    NB = int(os.environ.get('BOOT_NB', 20))


def closed_loops(alg, n):
    out = set()
    for p in nonbacktracking_paths(alg.D, n):
        if len(p) == n and not any(endpoint(p, alg.D)):
            c = alg.canonical(p)
            if len(c) == n:
                out.add(c)
    return sorted(out)


def loop_value(U, path):
    """Site-averaged (1/2) Tr of the product along path, lattice axes 0..D-1."""
    d = np.zeros(D, int)
    P = None
    for s in path:
        mu = abs(s) - 1
        if s < 0:
            d[mu] -= 1
        link = np.roll(U[mu], tuple(-d), axis=tuple(range(D)))
        if s < 0:
            link = qdag(link)
        else:
            d[mu] += 1
        P = link if P is None else qmul(P, link)
    return float(P[..., 0].mean())


def main():
    alg = Algebra(D)
    eqs = []
    for n in LENGTHS:
        for C in closed_loops(alg, n):
            eqs += alg.equations_for(C, BETA)
    uniq = {tuple(sorted((k, round(v, 12)) for k, v in e.items())): e for e in eqs}
    eqs = list(uniq.values())
    loops = sorted({k for e in eqs for k in e if k})
    nself = sum(1 for e in eqs for k in e if k and len(k) < 4)
    print(f'{len(eqs)} distinct equations from all loops of lengths {LENGTHS} (D={D}); '
          f'{len(loops)} distinct loops (max length {max(map(len, loops))})')
    lat = Lattice(D, L, BETA, seed=SEED)
    for _ in range(THERM):
        lat.sweep()
    vals = []
    for c in range(NCFG):
        for _ in range(SKIP):
            lat.sweep()
        vals.append([loop_value(lat.U, k) for k in loops])
    vals = np.array(vals)
    idx = {k: i for i, k in enumerate(loops)}
    A = np.zeros((len(eqs), len(loops)))
    b = np.zeros(len(eqs))
    for r, e in enumerate(eqs):
        for k, v in e.items():
            if k:
                A[r, idx[k]] += v
            else:
                b[r] += v
    res_cfg = vals @ A.T + b                                  # (cfg, eq)
    bins = np.array_split(res_cfg, NB)
    tot = sum(x.sum(0) for x in bins)
    full = tot / len(res_cfg)
    reps = np.array([(tot - x.sum(0)) / (len(res_cfg) - len(x)) for x in bins])
    err = np.sqrt((NB - 1) * np.mean((reps - reps.mean(0)) ** 2, 0))
    z = full / err
    scale = np.abs(A) @ np.abs(vals.mean(0)) + np.abs(b)
    print(f'residual / typical term size: median {np.median(np.abs(full) / scale):.2e}, '
          f'max {np.max(np.abs(full) / scale):.2e}')
    print(f'|z| > 3: {int(np.sum(np.abs(z) > 3))} of {len(z)} equations;  '
          f'max |z| = {np.max(np.abs(z)):.2f};  mean z^2 = {np.mean(z ** 2):.2f} '
          f'(about 1 if all equations hold; equations are correlated)')
    from scipy.stats import norm
    zN = norm.ppf(1 - 0.005 / len(z))
    okp = (np.max(np.abs(z)) < zN, np.mean(z ** 2) < 2)
    print(f"V1' (a) max|z| {np.max(np.abs(z)):.2f} < z_N {zN:.2f}: {okp[0]};  (b) mean z^2 < 2: {okp[1]}")
    worst = np.argsort(-np.abs(z))[:3]
    for r in worst:
        print(f'   worst eq {r}: residual {full[r]:+.2e} +- {err[r]:.1e}  (z = {z[r]:+.2f})')
    if D == 4 and len(eqs) > 5:
        print(f'   [V1 run-twice worst, eq 5]: z = {z[5]:+.2f}; terms: '
              + ', '.join(f'{v:+.3f}*w{k}' for k, v in sorted(eqs[5].items(), key=lambda kv: len(kv[0]))))
    ok = np.sum(np.abs(z) > 3) == 0         # as registered: EVERY equation
    print(f'V1 (every equation within 3 sigma, as registered) -> {"PASS" if ok else "FAIL"}')

    # negative control: a deliberately wrong equation must be caught
    bad = dict(eqs[0])
    k0 = next(k for k in bad if k)
    bad[k0] *= 1.02
    rb = vals @ np.array([bad.get(k, 0) for k in loops]) + bad.get((), 0)
    bb = np.array_split(rb, NB)
    tb = sum(x.sum() for x in bb)
    rp = np.array([(tb - x.sum()) / (len(rb) - len(x)) for x in bb])
    zb = rb.mean() / np.sqrt((NB - 1) * np.mean((rp - rp.mean()) ** 2))
    print(f'negative control (one coefficient off by 2%): z = {zb:+.1f}  '
          f'-> {"caught" if abs(zb) > 3 else "NOT caught (test has no power)"}')
    import json
    out_path = os.environ.get('BOOT_DUMP')
    if out_path:
        json.dump(dict(z=z.tolist(), res=full.tolist(), err=err.tolist(),
                       reps=reps.tolist(), control_z=float(zb), N=len(z)),
                  open(out_path, 'w'))
    print(f"V1' (c) control caught at |z| > 5: {abs(zb) > 5}  ->  V1' {'PASS' if okp[0] and okp[1] and abs(zb) > 5 else 'FAIL'}")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""POST-HOC (not registered), written after reading REPORT.txt.

1. K0 failed at kappa = 0, L = 10 with lambda(2) > lambda(1), the signature of
   largest-eigenvalue noise bias rather than a state. Check: (a) each operator's
   own normalised correlator (no maximisation, no bias); (b) the pure-noise
   floor of lambda_max built from the measured element errors (as SPEC-01).
2. The registered particle rule called kappa = 1 (4D, known glueball) "no
   particle" because L = 6 (L sqrt(sigma) ~ 2.3) shifts the mass, the effect
   ELIM-01 excluded L = 6 for. Re-fit with L >= 8 only.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'mirror'))
import analyse_mirror as AM  # noqa: E402

AM.LT = 16
rng = np.random.default_rng(3)


def load(L, k):
    z = np.load(os.path.join(HERE, 'raw', f'L{L}_k{k:.2f}.npz'))
    return [str(s) for s in z['names']], z['S1'], z['S2'], z['n']


print('1(a) single-operator C(t)/C(0) at kappa = 0 (unbiased); |z| > 3 flagged')
flags = 0
for L in (6, 8, 10):
    names, S1, S2, n = load(L, 0.0)
    for op in [s for s in names if s != 'poly']:
        g = AM.channel(names, S1, S2, n, [op])
        zs = [g['c'][t] / g['c_err'][t] for t in (1, 2)]
        flags += sum(abs(z) > 3 for z in zs)
        print(f"   L={L:2d} {op:6s}  C(1)/C(0) = {g['c'][1]:+.4f}({g['c_err'][1] * 1e4:.0f}) "
              f"z={zs[0]:+.1f}   C(2)/C(0) = {g['c'][2]:+.4f}({g['c_err'][2] * 1e4:.0f}) z={zs[1]:+.1f}")
print(f'   -> {flags} of 24 single-operator values beyond 3 sigma')

print('\n1(b) lambda_max vs its pure-noise floor (mean [95th pct]) at kappa = 0, 0.25, 0.5')
for k in (0.0, 0.25, 0.5):
    for L in (6, 8, 10):
        names, S1, S2, n = load(L, k)
        idx = [i for i, s in enumerate(names) if s != 'poly']
        C = AM.corr(S1, S2, n, idx)
        d, W = AM.gevp_basis(C)
        A = lambda C: np.einsum('ki,klt,lj->ijt', W, C * d[:, None, None] * d[None, :, None], W)
        full = A(C)
        reps = np.array([A(AM.corr(S1, S2, n, idx, drop=b)) for b in range(len(n))])
        sd = np.sqrt((len(n) - 1) * np.mean((reps - reps.mean(0)) ** 2, 0))
        cells = []
        for t in (1, 2):
            obs = np.linalg.eigvalsh(full[:, :, t])[-1]
            s = 0.5 * (sd[:, :, t] + sd[:, :, t].T)
            kk = s.shape[0]
            dr = []
            for _ in range(3000):
                X = rng.standard_normal((kk, kk)) * s
                X = np.triu(X) + np.triu(X, 1).T
                dr.append(np.linalg.eigvalsh(X)[-1])
            dr = np.array(dr)
            cells.append(f't={t}: {obs:.4f} vs noise {dr.mean():.4f} [{np.percentile(dr, 95):.4f}]'
                         + (' *' if obs > np.percentile(dr, 95) else '  '))
        print(f'   kappa={k:.2f} L={L:2d}  ' + '   '.join(cells))

print('\n2. particle rule with L >= 8 only (ELIM-01 convention), using m(0->1)')
for k in (0.75, 1.0):
    E, e = [], []
    for L in (8, 10):
        names, S1, S2, n = load(L, k)
        g = AM.channel(names, S1, S2, n, [s for s in names if s != 'poly'])
        E.append(g['m01'])
        e.append(g['m01_err'])
    E, e, Lv = np.array(E), np.array(e), np.array([8.0, 10.0])
    w = 1 / e ** 2
    c0 = np.sum(w * E) / np.sum(w)
    chi_c = np.sum(w * (E - c0) ** 2)
    kf = np.sum(w * E * Lv) / np.sum(w * Lv ** 2)
    ks = np.sum(w * E / Lv) / np.sum(w / Lv ** 2)
    print(f'   kappa={k:.2f}: m(0->1) = {E[0]:.3f}({e[0] * 1e3:.0f}), {E[1]:.3f}({e[1] * 1e3:.0f})'
          f'   const chi2 {chi_c:.1f}/1   flux k*L chi2 {np.sum(w * (E - kf * Lv) ** 2):.0f}/1'
          f'   scale-inv chi2 {np.sum(w * (E - ks / Lv) ** 2):.0f}/1')

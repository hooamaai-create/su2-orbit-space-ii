#!/usr/bin/env python3
"""POST-HOC diagnostic (not registered): is a flat lambda_max a state or noise?

Written after reading REPORT.txt, because in T1+, T1- and T2- the largest
principal correlator stops decaying at ~0.04, which makes m_eff(1->2) come out
lighter than the scalar. That is exactly what an accidental light state would
look like. It is also what the largest eigenvalue of a noisy K x K matrix
looks like: the maximum over many noisy directions is biased upward, by an
amount that grows with K, and jackknife errors cannot see a bias that every
sample shares.

Null model: a symmetric Gaussian matrix with the measured jackknife standard
deviation of every element of the whitened C(t), and zero mean. Its largest
eigenvalue is what lambda_max(t) reads when there is no signal at all.

Run from repo root:  python spectrum/noise_floor.py > spectrum/NOISE_FLOOR.txt
"""
import json
import os

import numpy as np

from analyse import CH, PRETTY, basis, corr, load

HERE = os.path.dirname(os.path.abspath(__file__))
NDRAW = 4000


def whitened(C, d, Q, W):
    Cn = C * d[:, None, None] * d[None, :, None]
    Cq = np.einsum('ak,abt,bl->klt', Q, Cn, Q)
    return np.einsum('ki,klt,lj->ijt', W, Cq, W)


def main():
    ch, L, _ = load()
    rng = np.random.default_rng(1)
    out = {}
    print('POST-HOC — observed lambda_max(t) vs pure-noise floor '
          '(mean [95th pct] of largest eigenvalue)')
    print(' chan  kept   t=2: observed  noise              '
          't=3: observed  noise              t=4: observed  noise')
    for c in CH:
        S1, S2, n = ch[c]
        C = corr(S1, S2, n)
        d, Q = basis(C)
        Cn = C * d[:, None, None] * d[None, :, None]
        Cq = np.einsum('ak,abt,bl->klt', Q, Cn, Q)
        w0, V0 = np.linalg.eigh(Cq[:, :, 0])
        W = V0 / np.sqrt(w0)
        A = whitened(C, d, Q, W)
        nb = len(n)
        reps = np.array([whitened(corr(S1, S2, n, drop=b), d, Q, W)
                         for b in range(nb)])
        sd = np.sqrt((nb - 1) * np.mean((reps - reps.mean(0)) ** 2, 0))
        k = A.shape[0]
        row = {}
        cells = []
        for t in (2, 3, 4):
            obs = float(np.linalg.eigvalsh(A[:, :, t])[-1])
            s = 0.5 * (sd[:, :, t] + sd[:, :, t].T)
            draws = []
            for _ in range(NDRAW):
                X = rng.standard_normal((k, k)) * s
                X = np.triu(X) + np.triu(X, 1).T
                draws.append(np.linalg.eigvalsh(X)[-1])
            draws = np.array(draws)
            mean, p95 = float(draws.mean()), float(np.percentile(draws, 95))
            row[t] = dict(observed=obs, noise_mean=mean, noise_p95=p95,
                          above_p95=bool(obs > p95))
            flag = ' *' if obs > p95 else '  '
            cells.append(f'{obs:7.4f}   {mean:.4f} [{p95:.4f}]{flag}')
        out[c] = dict(kept=k, t=row)
        print(f' {PRETTY[c]:4s}  {k:3d}    ' + '   '.join(cells))
    n_above = sum(v['above_p95'] for r in out.values() for v in r['t'].values())
    print(f'\n * = observed above the 95th percentile of pure noise. '
          f'{n_above} of {3 * len(CH)} cells (about {0.05 * 3 * len(CH):.1f} '
          f'expected by chance alone).')
    with open(os.path.join(HERE, 'noise_floor.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()

"""Semidefinite program: bounds on a Wilson loop from loop equations + positivity.

bounds(D, beta, K) maximises and minimises w(plaquette) over all assignments
of Wilson-loop values that
  - satisfy every loop equation whose loops all have length <= 2K, and
  - make every Gram block of paths of length <= K positive semidefinite.

Status of the numbers: they are what a floating-point interior-point solver
returns, accurate to the solver tolerance (~1e-7 here). They become *proven*
bounds only after a dual certificate is checked in exact arithmetic, which
this module does not do.
"""
import time

import cvxpy as cp
import numpy as np
import scipy.sparse as sp

from loops import PLAQ, Algebra, gram_blocks


def build(D, beta, K, verbose=False):
    t0 = time.time()
    alg = Algebra(D)
    blocks = gram_blocks(alg, K)
    lmax = 2 * K
    in_blocks = {k for b in blocks for row in b['keys'] for k in row}
    seeds = sorted(k for k in in_blocks if k and len(k) <= lmax - 4)
    eqs, seen = [], set()
    for C in seeds:
        for e in alg.equations_for(C, beta):
            if any(len(k) > lmax for k in e):
                continue
            key = tuple(sorted((k, round(v, 10)) for k, v in e.items()))
            if key not in seen:
                seen.add(key)
                eqs.append(e)
    loops = sorted({k for e in eqs for k in e if k} | {k for k in in_blocks if k})
    idx = {k: i + 1 for i, k in enumerate(loops)}         # 0 = the constant 1
    idx[()] = 0
    if verbose:
        sizes = sorted((len(b['paths']) for b in blocks), reverse=True)
        print(f'  D={D} K={K}: {len(blocks)} Gram blocks (largest {sizes[:4]}), '
              f'{len(eqs)} equations, {len(loops)} loop variables '
              f'[{time.time() - t0:.1f}s]')
    return alg, blocks, eqs, loops, idx


def bounds(D, beta, K, obs=PLAQ, verbose=False, solver='CLARABEL'):
    alg, blocks, eqs, loops, idx = build(D, beta, K, verbose)
    n = len(loops) + 1
    w = cp.Variable(n)
    cons = [w[0] == 1]
    rows, cols, vals = [], [], []
    for r, e in enumerate(eqs):
        for k, v in e.items():
            rows.append(r)
            cols.append(idx[k])
            vals.append(v)
    A = sp.csr_matrix((vals, (rows, cols)), shape=(len(eqs), n))
    cons.append(A @ w == 0)
    for b in blocks:
        keys = b['keys']
        m = len(keys)
        sel = sp.csr_matrix((np.ones(m * m), (np.arange(m * m),
                             [idx[keys[i][j]] for i in range(m) for j in range(m)])),
                            shape=(m * m, n))
        M = cp.reshape(sel @ w, (m, m), order='C')
        cons.append(0.5 * (M + M.T) >> 0)
    target = idx[alg.canonical(obs)]
    out = {}
    for sense in ('min', 'max'):
        prob = cp.Problem(cp.Minimize(w[target]) if sense == 'min'
                          else cp.Maximize(w[target]), cons)
        t0 = time.time()
        prob.solve(solver=solver)
        out[sense] = (float(w.value[target]) if w.value is not None else np.nan,
                      prob.status, round(time.time() - t0, 1))
    return out

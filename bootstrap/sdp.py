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


# ---------------------------------------------------------------------------
# Direct Clarabel formulation (no cvxpy): same program, ~100x less memory.
# PSD blocks enter in Clarabel's svec form: upper triangle, column-major,
# off-diagonals scaled by sqrt(2).
# ---------------------------------------------------------------------------

def _svec_index(m):
    rows, cols = np.triu_indices(m)
    order = np.lexsort((rows, cols))            # column-major over upper triangle
    rows, cols = rows[order], cols[order]
    scale = np.where(rows == cols, 1.0, np.sqrt(2.0))
    return rows, cols, scale


def bounds_direct(D, beta, K, obs=PLAQ, verbose=False, settings=None, built=None):
    import clarabel
    alg, blocks, eqs, loops, idx = built or build(D, beta, K, verbose)
    nvar = len(loops)                            # w_1..w_n; w_0 = 1 is a constant
    Arows, Acols, Avals, b, cones = [], [], [], [], []
    r0 = 0
    # equalities:  sum_k c_k w_k = -c_0   ->  A x + s = b, s in Zero cone
    for e in eqs:
        for k, v in e.items():
            if k:
                Arows.append(r0)
                Acols.append(idx[k] - 1)
                Avals.append(v)
        b.append(-e.get((), 0.0))
        r0 += 1
    cones.append(clarabel.ZeroConeT(len(eqs)))
    # PSD blocks:  svec(M(w)) = s,  A x + s = b  with A = -svec(E_k), b = svec(M_const)
    for blk in blocks:
        keys = blk['keys']
        m = len(keys)
        rr, cc, sc = _svec_index(m)
        for t, (i, j, f) in enumerate(zip(rr, cc, sc)):
            k = keys[i][j]
            if k:
                Arows.append(r0 + t)
                Acols.append(idx[k] - 1)
                Avals.append(-f)
                b.append(0.0)
            else:
                b.append(f)
        r0 += len(rr)
        cones.append(clarabel.PSDTriangleConeT(m))
    A = sp.csc_matrix((Avals, (Arows, Acols)), shape=(r0, nvar))
    b = np.array(b)
    P = sp.csc_matrix((nvar, nvar))
    tgt = idx[alg.canonical(obs)] - 1
    out = {}
    for sense, sgn in (('min', 1.0), ('max', -1.0)):
        q = np.zeros(nvar)
        q[tgt] = sgn
        st = clarabel.DefaultSettings()
        st.verbose = False
        for kk, vv in (settings or {}).items():
            setattr(st, kk, vv)
        t0 = time.time()
        sol = clarabel.DefaultSolver(P, q, A, b, cones, st).solve()
        val = float(sol.x[tgt]) if sol.x is not None else np.nan
        out[sense] = (val, str(sol.status), round(time.time() - t0, 1))
    return out


def bounds_scs(D, beta, K, obs=PLAQ, verbose=False, eps=1e-7, max_iters=200000,
               built=None):
    """Same program via SCS (first-order, low memory).

    SCS vectorises a PSD block as its lower triangle, column-major, with
    off-diagonals scaled by sqrt(2); for a symmetric block that is the upper
    triangle row-major, which is what is built here.
    """
    import scs
    alg, blocks, eqs, loops, idx = built or build(D, beta, K, verbose)
    nvar = len(loops)
    Ar, Ac, Av, b = [], [], [], []
    r0 = 0
    for e in eqs:
        for k, v in e.items():
            if k:
                Ar.append(r0)
                Ac.append(idx[k] - 1)
                Av.append(v)
        b.append(-e.get((), 0.0))
        r0 += 1
    nz = r0
    sizes = []
    for blk in blocks:
        keys = blk['keys']
        m = len(keys)
        t = 0
        for j in range(m):                 # lower triangle, column-major
            for i in range(j, m):
                f = 1.0 if i == j else np.sqrt(2.0)
                k = keys[i][j]
                if k:
                    Ar.append(r0 + t)
                    Ac.append(idx[k] - 1)
                    Av.append(-f)
                    b.append(0.0)
                else:
                    b.append(f)
                t += 1
        r0 += t
        sizes.append(m)
    A = sp.csc_matrix((Av, (Ar, Ac)), shape=(r0, nvar))
    b = np.array(b)
    tgt = idx[alg.canonical(obs)] - 1
    out = {}
    for sense, sgn in (('min', 1.0), ('max', -1.0)):
        c = np.zeros(nvar)
        c[tgt] = sgn
        t0 = time.time()
        solver = scs.SCS(dict(A=A, b=b, c=c), dict(z=nz, s=sizes),
                         eps_abs=eps, eps_rel=eps, max_iters=max_iters, verbose=False)
        sol = solver.solve()
        out[sense] = (float(sol['x'][tgt]), sol['info']['status'],
                      round(time.time() - t0, 1))
    return out

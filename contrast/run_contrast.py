#!/usr/bin/env python3
"""Run the 2D-vs-4D contrast and write contrast/results.json.

Same engine, same observables, both dimensions:
  plaquette, Wilson loops W(R,T), Creutz ratios, slice correlator C(t).
Every 2D number has an exact answer to fail against; every 4D number is set
against leading strong coupling so the point where 2D-style reasoning stops
working is measured, not asserted.

Run from repo root:  python contrast/run_contrast.py
"""
import json
import time

import numpy as np

from su2_lattice import Lattice, u_exact

RMAX = 4
RUNS = (
    [dict(D=2, L=32, beta=b, therm=200, meas=1200) for b in (1.0, 2.0, 4.0, 8.0)]
    + [dict(D=4, L=8, beta=b, therm=150, meas=700)
       for b in (1.0, 1.5, 2.0, 2.2, 2.4, 2.6)]
)
WILSON_EVERY = 4
NBINS = 20


def jack(samples, f):
    """Binned jackknife of f(mean over samples along axis 0)."""
    x = np.asarray(samples)
    bins = np.array_split(x, NBINS)
    tot = sum(b.sum(0) for b in bins)
    n = len(x)
    full = f(tot / n)
    reps = np.array([f((tot - b.sum(0)) / (n - len(b))) for b in bins])
    err = np.sqrt((NBINS - 1) * np.mean((reps - reps.mean(0)) ** 2, axis=0))
    return full, err


def creutz(W):
    """chi(R) = -ln W(R,R) W(R-1,R-1) / W(R,R-1)^2 for R = 2..RMAX."""
    out = []
    for r in range(2, W.shape[0] + 1):
        num = W[r - 1, r - 1] * W[r - 2, r - 2]
        den = W[r - 1, r - 2] * W[r - 2, r - 1]
        out.append(-np.log(num / den) if num > 0 and den > 0 else np.nan)
    return np.array(out)


def slice_corr(ops):
    """C(t) = <O(t0) O(t0+t)> - <O>^2 averaged over t0; returns components."""
    ops = np.asarray(ops)                    # (n, L)
    L = ops.shape[1]
    prod = np.stack([np.mean(ops * np.roll(ops, -t, axis=1), axis=1)
                     for t in range(L // 2 + 1)], axis=1)
    mean = ops.mean(1)
    return np.concatenate([prod, mean[:, None]], axis=1)


def corr_from(v):
    L2 = len(v) - 1
    return v[:L2] - v[L2] ** 2


def main():
    out = []
    t_all = time.time()
    for cfg in RUNS:
        t0 = time.time()
        D, L, beta = cfg['D'], cfg['L'], cfg['beta']
        lat = Lattice(D, L, beta, seed=int(100 * beta) + D)
        for _ in range(cfg['therm']):
            lat.sweep()
        P, Wl, S = [], [], []
        for i in range(cfg['meas']):
            lat.sweep()
            P.append(lat.plaquette())
            S.append(lat.slice_op())
            if i % WILSON_EVERY == 0:
                Wl.append(lat.wilson(RMAX).ravel())
        u = u_exact(beta)

        p, pe = jack(P, lambda m: m)
        W, We = jack(Wl, lambda m: m.reshape(RMAX, RMAX))
        chi, chie = jack(Wl, lambda m: creutz(m.reshape(RMAX, RMAX)))
        comp = slice_corr(S)
        C, Ce = jack(comp, corr_from)
        meff, meffe = jack(comp, lambda m: np.log(np.abs(corr_from(m)[0] /
                                                           corr_from(m)[1])))
        rec = dict(
            D=D, L=L, beta=beta, therm=cfg['therm'], meas=cfg['meas'],
            u=u,
            plaquette=[float(p), float(pe)],
            W=W.tolist(), W_err=We.tolist(),
            W_exact2D=[[u ** (r * t) for t in range(1, RMAX + 1)]
                       for r in range(1, RMAX + 1)],
            creutz=chi.tolist(), creutz_err=chie.tolist(),
            sigma_strong=-np.log(u),
            C=C.tolist(), C_err=Ce.tolist(),
            C_rel=(C / C[0]).tolist(), C_rel_err=(Ce / abs(C[0])).tolist(),
            meff01=[float(meff), float(meffe)],
            m_strong=-4 * np.log(u),
            seconds=round(time.time() - t0, 1),
        )
        out.append(rec)
        print(f"D={D} L={L} beta={beta:3.1f}  <P>={p:.4f}({pe * 1e4:.0f}) "
              f"u={u:.4f}  chi2={chi[0]:.3f}({chie[0] * 1e3:.0f}) "
              f"-ln u={-np.log(u):.3f}  C(1)/C(0)={C[1] / C[0]:+.4f}"
              f"({Ce[1] / abs(C[0]) * 1e4:.0f})  "
              f"m01={meff:.2f}({meffe * 100:.0f})  -4ln u={-4 * np.log(u):.2f}"
              f"  [{rec['seconds']}s]", flush=True)
    with open('contrast/results.json', 'w') as f:
        json.dump(dict(runs=out, rmax=RMAX, nbins=NBINS,
                       wall_seconds=round(time.time() - t_all, 1)), f, indent=1)
    print("wrote contrast/results.json")


if __name__ == '__main__':
    main()

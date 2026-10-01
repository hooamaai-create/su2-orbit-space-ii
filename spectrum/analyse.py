#!/usr/bin/env python3
"""SPEC-01 analysis: score S1-S4 from spectrum/raw/stream*.npz.

Written and committed before any stream finished. Choices, all from SPEC.md:
  - 20 jackknife bins (5 per stream x 4 streams);
  - per channel: normalise by diag C(0), prune modes of C(0) below 1e-3 of the
    largest (basis fixed from the full ensemble, reused in every jackknife
    sample), whiten at t0 = 0, principal correlators = eigenvalues of the
    whitened C(t) at each t;
  - ground state = largest principal correlator; cosh effective mass;
    primary A1+ mass = m_eff(1->2), as in ELIM-01;
  - light-state bound w_max(E) = min_{t>=1} [lambda_max(t) + 2 sigma] / f_E(t),
    f_E(t) = cosh(E(t - L/2)) / cosh(E L/2).

Run from repo root:  python spectrum/analyse.py > spectrum/REPORT.txt
"""
import glob
import json
import os

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
CH = ['A1p', 'A1m', 'A2p', 'A2m', 'Ep', 'Em', 'T1p', 'T1m', 'T2p', 'T2m']
PRETTY = {c: c[:-1] + ('+' if c[-1] == 'p' else '-') for c in CH}
CONTINUUM = dict(A1p='0++', A1m='0-+', A2p='3++', A2m='3-+', Ep='2++',
                 Em='2-+', T1p='1++', T1m='1-+', T2p='2++', T2m='2-+')
ELIM_A1P = (1.197, 0.117)          # ELIM-01, beta 2.3, L 10, m_eff(1->2)
PRUNE = 1e-3
bar = '=' * 78


def load():
    files = sorted(glob.glob(os.path.join(HERE, 'raw', 'stream*.npz')))
    data = [dict(np.load(f)) for f in files]
    L = int(data[0]['L'])
    ch = {}
    for c in CH:
        S1 = np.concatenate([d[f'{c}_S1'] for d in data])
        S2 = np.concatenate([d[f'{c}_S2'] for d in data])
        n = np.concatenate([d[f'{c}_n'] for d in data])
        ch[c] = (S1, S2, n)
    return ch, L, len(files)


def corr(S1, S2, n, drop=None):
    keep = np.ones(len(n), bool)
    if drop is not None:
        keep[drop] = False
    N = n[keep].sum()
    Y = S1[keep].sum(0) / N
    M = S2[keep].sum(0) / N
    C = M - np.einsum('a,b->ab', Y, Y)[:, :, None]
    return 0.5 * (C + C.transpose(1, 0, 2))


def basis(C):
    d = 1.0 / np.sqrt(np.diag(C[:, :, 0]))
    Cn = C * d[:, None, None] * d[None, :, None]
    w, V = np.linalg.eigh(Cn[:, :, 0])
    keep = w > PRUNE * w.max()
    return d, V[:, keep]


def principal(C, d, Q):
    Cn = C * d[:, None, None] * d[None, :, None]
    Cq = np.einsum('ak,abt,bl->klt', Q, Cn, Q)
    w0, V0 = np.linalg.eigh(Cq[:, :, 0])
    W = V0 / np.sqrt(np.maximum(w0, 1e-300))
    lam = np.array([np.linalg.eigvalsh(W.T @ Cq[:, :, t] @ W)[::-1]
                    for t in range(Cq.shape[2])])
    return lam                              # (T, k), descending


def meff(c, L):
    out = []
    h = L / 2
    for t in range(len(c) - 1):
        if c[t] <= 0 or c[t + 1] <= 0 or c[t] <= c[t + 1] or h - t - 1 <= 0:
            out.append(np.nan)
            continue
        r = c[t] / c[t + 1]
        try:
            out.append(brentq(lambda m: np.cosh(m * (h - t)) /
                              np.cosh(m * (h - t - 1)) - r, 1e-6, 30.0))
        except ValueError:
            out.append(np.nan)
    return np.array(out)


def f_E(E, t, L):
    return np.cosh(E * (t - L / 2)) / np.cosh(E * L / 2)


def analyse_channel(S1, S2, n, L):
    C = corr(S1, S2, n)
    d, Q = basis(C)
    lam = principal(C, d, Q)[:, 0]
    m = meff(lam, L)
    nb = len(n)
    reps_l, reps_m = [], []
    for b in range(nb):
        lb = principal(corr(S1, S2, n, drop=b), d, Q)[:, 0]
        reps_l.append(lb)
        reps_m.append(meff(lb, L))
    reps_l, reps_m = np.array(reps_l), np.array(reps_m)
    jk = lambda r: np.sqrt((nb - 1) * np.nanmean((r - np.nanmean(r, 0)) ** 2, 0))
    return dict(K=int(C.shape[0]), kept=int(Q.shape[1]),
                lam=lam.tolist(), lam_err=jk(reps_l).tolist(),
                meff=m.tolist(), meff_err=jk(reps_m).tolist())


def w_max(r, E, L):
    lam, err = np.array(r['lam']), np.array(r['lam_err'])
    t = np.arange(1, len(lam))
    vals = (lam[1:] + 2 * err[1:]) / f_E(E, t, L)
    k = int(np.argmin(vals))
    return float(min(vals[k], 1.0)), int(t[k])


def fmt(x, e):
    if not np.isfinite(x):
        return '   n/a   '
    return f'{x:.3f}({e * 1e3:.0f})' if np.isfinite(e) else f'{x:.3f}(?)'


def main():
    ch, L, ns = load()
    R = {c: analyse_channel(*ch[c], L) for c in CH}
    m0, m0e = R['A1p']['meff'][1], R['A1p']['meff_err'][1]
    n_cfg = int(ch['A1p'][2].sum())

    print(bar)
    print(f'SPEC-01  beta 2.3, L = {L}, {ns} streams, {n_cfg} measurements')
    print(bar)
    print(' channel  (J^PC)  ops kept   m_eff(0->1)   m_eff(1->2)   '
          'lambda_max(t = 1, 2, 3)')
    for c in CH:
        r = R[c]
        lam = ', '.join(f"{l:.4f}({e * 1e4:.0f})" for l, e in
                        zip(r['lam'][1:4], r['lam_err'][1:4]))
        print(f" {PRETTY[c]:4s}    ({CONTINUUM[c]})  {r['K']:3d} {r['kept']:3d}   "
              f"{fmt(r['meff'][0], r['meff_err'][0]):12s}  "
              f"{fmt(r['meff'][1], r['meff_err'][1]):12s}  {lam}")

    out = dict(channels=R, L=L, streams=ns, n_cfg=n_cfg)

    print('\n' + bar)
    print('S1 — control')
    print(bar)
    dz = abs(m0 - ELIM_A1P[0]) / np.hypot(m0e, ELIM_A1P[1])
    s1a = dz < 2
    wA, tA = w_max(R['A1p'], m0, L)
    s1b = wA >= 0.5
    print(f"  (a) A1+ m_eff(1->2) = {m0:.3f}({m0e * 1e3:.0f})  vs ELIM-01 "
          f"{ELIM_A1P[0]:.3f}({ELIM_A1P[1] * 1e3:.0f}): {dz:.2f} sigma -> "
          f"{'PASS' if s1a else 'FAIL'}")
    print(f"  (b) bound must NOT exclude the 0++ that exists: "
          f"w_max(E = m0) in A1+ = {wA:.3f} (at t={tA}) -> "
          f"{'PASS' if s1b else 'FAIL (method broken)'}")
    s1 = s1a and s1b
    out['S1'] = dict(a=bool(s1a), dz=float(dz), b=bool(s1b), w=wA,
                     pass_=bool(s1))

    print('\n' + bar)
    print('S2 — rotational check: E+ vs T2+ (both 2++)')
    print(bar)
    e, ee = R['Ep']['meff'][0], R['Ep']['meff_err'][0]
    t2, t2e = R['T2p']['meff'][0], R['T2p']['meff_err'][0]
    z2 = abs(e - t2) / np.hypot(ee, t2e)
    rel = abs(e - t2) / np.mean([e, t2])
    s2 = z2 < 2 or rel < 0.15
    print(f"  E+ {fmt(e, ee)}   T2+ {fmt(t2, t2e)}   {z2:.1f} sigma, "
          f"{rel:.1%} apart -> {'PASS' if s2 else 'WARNING: lattice artefact'}")
    out['S2'] = dict(z=float(z2), rel=float(rel), pass_=bool(s2))

    print('\n' + bar)
    print('S3 — is any channel lighter than the scalar?  (m_eff(0->1))')
    print(bar)
    a0, a0e = R['A1p']['meff'][0], R['A1p']['meff_err'][0]
    s3 = []
    for c in CH[1:]:
        m, me = R[c]['meff'][0], R[c]['meff_err'][0]
        ok = bool(np.isfinite(m) and m > a0 - 2 * np.hypot(a0e, me))
        s3.append(ok)
        print(f"  {PRETTY[c]:4s} {fmt(m, me):12s}  vs A1+ {fmt(a0, a0e)}  "
              f"ratio {m / a0:5.2f}  -> {'heavier' if ok else 'NOT heavier'}")
    s3p = all(s3)
    print(f"  S3 -> {'PASS: the scalar is the lightest channel' if s3p else 'FAIL'}"
          f"  ({sum(s3)}/{len(s3)})")
    out['S3'] = dict(each=s3, pass_=bool(s3p))

    print('\n' + bar)
    print('S4 — the loophole bound: max weight a light state can have in the basis')
    print(bar)
    print(f"  thresholds: massless (E=0) < 0.10, sub-half-scalar "
          f"(E = m0/2 = {m0 / 2:.3f}) < 0.30")
    s4 = {}
    for c in CH:
        w0, t0 = w_max(R[c], 0.0, L)
        wh, th = w_max(R[c], m0 / 2, L)
        ok = w0 < 0.10 and wh < 0.30
        s4[c] = dict(w0=w0, t0=t0, whalf=wh, thalf=th, closed=bool(ok))
        print(f"  {PRETTY[c]:4s}  w_max(0) = {w0:.3f} (t={t0})   "
              f"w_max(m0/2) = {wh:.3f} (t={th})   -> "
              f"{'closed' if ok else 'OPEN'}")
    n_closed = sum(v['closed'] for v in s4.values())
    verdict = ('LOOPHOLE CLOSED in all 10 channels (basis-visible states)'
               if n_closed == 10 else
               f'LOOPHOLE OPEN in {10 - n_closed} channel(s): '
               + ', '.join(PRETTY[c] for c, v in s4.items() if not v['closed']))
    if not s1:
        verdict = 'NOT SCORED — S1 control failed (instrument)'
    print(f"  S4 -> {verdict}")
    out['S4'] = dict(channels=s4, n_closed=int(n_closed), verdict=verdict)

    with open(os.path.join(HERE, 'results.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()

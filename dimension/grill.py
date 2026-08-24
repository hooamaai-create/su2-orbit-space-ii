#!/usr/bin/env python3
"""Adversarial self-review of the DIM-01 audit.

DIM-01 attacked the model. This attacks DIM-01. Five objections, each with a
run behind it rather than an opinion:

  G1  Is the ensemble "dispersion" real, or is it sampling noise in the
      per-probe estimator? Split-half reliability answers this.
  G2  Does the dispersion do any work at all? Re-run with the spread deleted
      (every probe set to the ensemble mean). If the flow still passes through
      the upper threshold, the "dispersion carries it past 4" claim -- the one
      piece of the model DIM-01 said survived -- is wrong.
  G3  Does the softening width w carry it instead? Push w toward a hard step
      and see whether the flow stalls at 4. DIM-01 swept w over 0.25-1.00 and
      never approached zero, which is where the interesting behaviour is.
  G4  Is the halting value pinned by the rate, or by the stopping tolerance?
      Vary gamma_stop over four decades.
  G5  Is the halt a property of the geometry or of the measuring window? Vary
      the window the dynamics measures in, and inspect the frozen states.

Run:  PYTHONPATH=dimension python dimension/grill.py --out dimension/results/grill.json
"""
import argparse
import json
import time

import numpy as np

import dim_flow
from dim_flow import (AMBIENT, apply_cal, calibration, contract, dim_vs_scale,
                      init_state, instability_rate, ir_extent, measure,
                      resample_like, run_flow)

BAR = '=' * 70


# --------------------------------------------------------------------------
# G1 — how much of the per-probe spread is signal?
# --------------------------------------------------------------------------


def split_half(x, rng, nprobe=600, f_r=dim_flow.F_R, minpts=16, chunk=100):
    """Estimate each probe's dimension twice from disjoint halves of its own
    neighbours. The correlation between the two halves is the reliability of a
    half-length estimate; Spearman-Brown gives the full-length reliability.

    Reliability is the fraction of the observed variance that is signal. If it
    is near zero, the "ensemble dispersion" is the estimator's own noise and
    every statement built on the spread is a statement about sampling error.
    """
    n = len(x)
    r = f_r * ir_extent(x)
    idx = rng.choice(n, size=min(nprobe, n), replace=False)
    a_half, b_half, full, counts = [], [], [], []
    for a in range(0, len(idx), chunk):
        sub = x[idx[a:a + chunk]]
        dd = np.linalg.norm(sub[:, None, :] - x[None, :, :], axis=2)
        for row in dd:
            t = row[(row > 0) & (row < r)]
            if len(t) < minpts:
                continue
            lg = np.log(r / t)
            perm = rng.permutation(len(t))
            h1, h2 = lg[perm[:len(t) // 2]], lg[perm[len(t) // 2:]]
            if h1.sum() <= 0 or h2.sum() <= 0:
                continue
            a_half.append(len(h1) / h1.sum())
            b_half.append(len(h2) / h2.sum())
            full.append(len(lg) / lg.sum())
            counts.append(len(t))
    a_half, b_half = np.array(a_half), np.array(b_half)
    full, counts = np.array(full), np.array(counts)
    if len(full) < 20:
        return None
    r_half = float(np.corrcoef(a_half, b_half)[0, 1])
    rel = 2 * r_half / (1 + r_half) if r_half > -1 else 0.0
    rel = float(np.clip(rel, 0.0, 1.0))
    return dict(n_probes=len(full), median_neighbours=float(np.median(counts)),
                observed_sd=float(np.std(full)),
                signal_sd=float(np.std(full) * np.sqrt(rel)),
                noise_sd=float(np.std(full) * np.sqrt(1 - rel)),
                half_corr=r_half, reliability=rel,
                predicted_noise_sd=float(np.mean(full) /
                                         np.sqrt(np.mean(counts))))


def u_of(dloc, w, d_collapse=4.0):
    return float(np.mean(1.0 / (1.0 + np.exp(-(dloc - d_collapse) / w))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='dimension/results/grill.json')
    ap.add_argument('--n', type=int, default=3000)
    args = ap.parse_args()

    t0 = time.time()
    cal = calibration(args.n)
    out = {}

    # ---------------------------------------------------------------- G1/G2
    print(BAR)
    print('G1 — is the "ensemble dispersion" signal or estimator noise?')
    print(BAR)
    r, x = run_flow(cal, seed=0, n=args.n, max_steps=1200,
                    capture_at=(6.0, 4.5, 3.2))
    print(f"  reference flow: {r['d_initial']:.3f} -> {r['d_final']:.3f} in "
          f"{r['steps']} steps ({r['halt']})")
    g1 = []
    for target, state in sorted(r['captures'].items(), reverse=True):
        rng = np.random.default_rng(11)
        sh = split_half(state, rng)
        if sh is None:
            continue
        d_here = apply_cal(measure(state, np.random.default_rng(12), 400)[0],
                           cal)
        sh.update(target=target, d_measured=float(d_here))
        g1.append(sh)
        print(f"  at d = {d_here:5.3f}: per-probe SD {sh['observed_sd']:.3f} "
              f"= signal {sh['signal_sd']:.3f} + noise {sh['noise_sd']:.3f}"
              f"   (reliability {sh['reliability']:.3f}, "
              f"{sh['median_neighbours']:.0f} neighbours/probe)")
    out['G1_split_half'] = g1

    print()
    print(BAR)
    print('G2 — does the dispersion do any work? Delete it and re-run.')
    print(BAR)
    # what fraction of U survives when every probe is set to the mean?
    st = r['captures'].get(3.2)
    if st is not None:
        _, dl_raw, _ = measure(st, np.random.default_rng(13), 800)
        dl = apply_cal(dl_raw, cal)
        u_real = u_of(dl, 0.50)
        u_flat = u_of(np.full_like(dl, dl.mean()), 0.50)
        print(f"  near the floor (mean local d = {dl.mean():.3f}):")
        print(f"    U with the real spread      {u_real:.4f}")
        print(f"    U with the spread deleted   {u_flat:.4f}"
              f"   ({100 * u_flat / u_real:.0f}% of it survives)")
        out['G2_U_decomposition'] = dict(mean_local_d=float(dl.mean()),
                                         u_with_spread=u_real,
                                         u_flat=u_flat)
    rf, _ = run_flow(cal, seed=0, n=args.n, max_steps=1200,
                     no_dispersion=True)
    print(f"  flow with dispersion deleted: {rf['d_initial']:.3f} -> "
          f"{rf['d_final']:.3f} in {rf['steps']} steps ({rf['halt']})")
    print(f"    -> passes through the upper threshold: "
          f"{'YES' if rf['d_final'] < 3.6 else 'NO, stalls'}")
    out['G2_no_dispersion'] = dict(d_initial=rf['d_initial'],
                                   d_final=rf['d_final'], steps=rf['steps'],
                                   halt=rf['halt'])

    # ------------------------------------------------------------------- G3
    print()
    print(BAR)
    print('G3 — is it the softening width w that carries the flow past 4?')
    print(BAR)
    g3 = []
    for w in (0.500, 0.250, 0.125, 0.060, 0.030):
        rw, _ = run_flow(cal, seed=0, n=args.n, w=w, max_steps=1200)
        g3.append(dict(w=w, d_final=rw['d_final'], steps=rw['steps'],
                       halt=rw['halt']))
        print(f"  w = {w:.3f}   halt {rw['d_final']:6.3f} after "
              f"{rw['steps']:4d} steps ({rw['halt']})")
    out['G3_hard_threshold'] = g3

    print()
    print(BAR)
    print('G3b — the crossed control: remove BOTH the spread and the softening')
    print(BAR)
    g3b = []
    for w, nd in ((0.030, False), (0.500, True), (0.030, True), (0.010, True)):
        rc, _ = run_flow(cal, seed=0, n=args.n, w=w, no_dispersion=nd,
                         max_steps=1200)
        g3b.append(dict(w=w, no_dispersion=nd, d_final=rc['d_final'],
                        steps=rc['steps'], halt=rc['halt']))
        print(f"  w = {w:.3f}  dispersion "
              f"{'deleted' if nd else 'intact ':8s} -> halt "
              f"{rc['d_final']:6.3f} in {rc['steps']:4d} steps ({rc['halt']})")
    out['G3b_crossed'] = g3b

    # ------------------------------------------------------------------- G4
    print()
    print(BAR)
    print('G4 — is the halt pinned by the rate or by the stopping tolerance?')
    print(BAR)
    g4 = []
    for gs in (1e-2, 1e-3, 1e-4, 1e-5):
        rg, _ = run_flow(cal, seed=0, n=args.n, gamma_stop=gs, max_steps=3000)
        g4.append(dict(gamma_stop=gs, d_final=rg['d_final'],
                       steps=rg['steps'], halt=rg['halt']))
        print(f"  gamma_stop = {gs:.0e}   halt {rg['d_final']:6.3f} after "
              f"{rg['steps']:4d} steps ({rg['halt']})")
    out['G4_tolerance'] = g4

    # ------------------------------------------------------------------- G5
    print()
    print(BAR)
    print('G5 — is the halt a property of the geometry or of the window?')
    print(BAR)
    g5 = []
    for f_r in (1.0, 1.4, 2.0):
        cal_f = calibration(args.n, f_r=f_r)
        rr, xx = run_flow(cal_f, seed=0, n=args.n, f_r=f_r, max_steps=1200)
        rng = np.random.default_rng(97)
        big = resample_like(xx, 20000, rng)
        prof = dim_vs_scale(big, cal, rng, fracs=(0.5, 1.0, 1.4, 2.0))
        line = '  '.join(f"{p['f_r']:.1f}:{p['d']:5.2f}" for p in prof)
        g5.append(dict(f_r=f_r, d_final=rr['d_final'], steps=rr['steps'],
                       halt=rr['halt'], frozen_profile=prof))
        print(f"  dynamics measures at R = {f_r:.1f}: halts at "
              f"{rr['d_final']:.3f}  |  that frozen state read at "
              f"R = {line}")
    out['G5_window'] = g5

    out['wall_seconds'] = round(time.time() - t0, 1)
    with open(args.out, 'w') as f:
        json.dump(out, f, indent=1)
    print(f"\nwrote {args.out}  ({out['wall_seconds']} s)")


if __name__ == '__main__':
    main()

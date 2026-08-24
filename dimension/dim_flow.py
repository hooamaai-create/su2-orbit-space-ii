#!/usr/bin/env python3
"""Dimension as a MEASURED statistic of correlations: the decay flow, and where it halts.

Independent re-implementation of the toy model, written to test one claim:

    "nothing in the code said 'go to 3'; the number emerged."

Nothing here stores a dimension. The state is a cloud of N correlated events in
an ambient space; the effective dimension is *read off* it at every step by a
fixed-radius maximum-likelihood estimator, exactly as one would read it off
data. The estimator is free to return anything in [0, ambient] and the loop is
free to run to zero.

What IS in the code, stated plainly so the emergence claim can be audited:

  d_collapse   the bound-state / circular-orbit stability threshold. In d
               spatial dimensions the two-body effective potential is
               -k/r^(d-2) + L^2/(2 mu r^2); circular orbits are stable iff
               d < 4, and the Coulomb problem falls to the centre for d >= 4.
               This is the (d-4) analysis. It enters as the point above which
               structure is erased -- i.e. it sets where the decay STARTS.
  d_grav       the dimension below which the long-range interaction that does
               the erasing carries no propagating degrees of freedom. It
               enters as the zero of the rate's coefficient -- i.e. it sets
               where the decay can STOP.

Both are inputs. The emergence claim is therefore not "3 was never mentioned";
it is the weaker, checkable claim that the flow halts AT THE LOWER THRESHOLD
rather than stalling at the upper one, overshooting it, or drifting with the
soft parameter w. dim_sweep.py tests exactly that, including the control that
matters most: move d_grav and see whether the halt follows it.

Run:  python dimension/dim_flow.py --out dimension/results/flow_default.json
"""
import argparse
import json
import time

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

AMBIENT = 12
F_R = 1.4          # measurement ball radius, in units of the cloud's IR extent
THETA = 0.15       # a principal axis is "resolvable" above this x IR extent
MINPTS = 8         # probes with fewer neighbours in the ball are discarded


# --------------------------------------------------------------------------
# state
# --------------------------------------------------------------------------


def init_state(rng, n, ambient, decay):
    """Anisotropic Gaussian cloud; principal axis k has std decay**k.

    `decay` fixes how many directions are resolvable inside the measurement
    ball and therefore the INITIAL measured dimension. It is a property of the
    initial condition; the dynamics never reads it.
    """
    return rng.standard_normal((n, ambient)) * decay ** np.arange(ambient)


def ir_extent(x):
    """The cloud's largest principal std: the scale the ruler is tied to.

    Tying the measurement ball to the IR extent rather than to a quantile of
    the pair distribution is what keeps the ruler from collapsing along with
    the microstructure it is measuring. With a quantile-fixed window the flow
    has a spurious feedback (dying axes concentrate pairs, the window shrinks,
    the surviving axes look isotropic again) and settles into a steady state
    near d = 4.5 that is an artifact of the estimator, not of the physics.
    """
    c = x - x.mean(0)
    return float(np.sqrt(np.linalg.eigvalsh(c.T @ c / len(c)).max()))


# --------------------------------------------------------------------------
# the estimator: dimension is measured, never stored
# --------------------------------------------------------------------------


def measure(x, rng, nprobe=400, f_r=F_R, minpts=MINPTS, chunk=100):
    """Fixed-radius (Hill / Levina-Bickel) MLE of the intrinsic dimension.

    For a point process on a d-dimensional set, the neighbour distances T_j
    inside a ball of radius R satisfy d_hat = m / sum_j log(R/T_j). Returns the
    pooled global estimate and the per-probe values, whose spread IS the
    ensemble dispersion the rate needs -- the flow is driven by the fraction of
    the ensemble that is locally super-critical, not by the mean.
    """
    n = len(x)
    r = f_r * ir_extent(x)
    idx = rng.choice(n, size=min(nprobe, n), replace=False)
    num = den = 0.0
    per = []
    for a in range(0, len(idx), chunk):
        sub = x[idx[a:a + chunk]]
        dd = np.linalg.norm(sub[:, None, :] - x[None, :, :], axis=2)
        for row in dd:
            t = row[(row > 0) & (row < r)]
            if len(t) < minpts:
                continue
            s = float(np.log(r / t).sum())
            if s <= 0:
                continue
            num += len(t)
            den += s
            per.append(len(t) / s)
    if den <= 0 or not per:
        return float('nan'), np.array([]), r
    return num / den, np.asarray(per), r


def calibration(n, ambient=AMBIENT, f_r=F_R, nprobe=400, kmax=8, seed=7):
    """Instrument calibration: measured dimension vs known dimension.

    The MLE has a mild multiplicative bias on a Gaussian (the density is not
    uniform inside the ball). It is measured here on isotropic clouds of known
    dimension at the same N and window as the flow, fitted linearly, and the
    inverse is applied to every dimension this module reports. Without it, a
    halt "at 2.97" would be uninterpretable: a genuinely 3-dimensional cloud
    reads about 2.66 raw.
    """
    ks, ds = [], []
    for k in range(1, kmax + 1):
        rng = np.random.default_rng(seed)
        x = np.zeros((n, ambient))
        x[:, :k] = rng.standard_normal((n, k))
        d, _, _ = measure(x, np.random.default_rng(seed + 1), nprobe, f_r)
        if np.isfinite(d):
            ks.append(float(k))
            ds.append(float(d))
    a, b = np.polyfit(ks, ds, 1)          # d_measured = a*k + b
    resid = float(np.max(np.abs(np.array(ds) - (a * np.array(ks) + b))))
    return dict(slope=float(a), intercept=float(b), max_resid=resid,
                k=ks, d_raw=ds)


def apply_cal(d_raw, cal):
    return (d_raw - cal['intercept']) / cal['slope']


# --------------------------------------------------------------------------
# the rate
# --------------------------------------------------------------------------


def instability_rate(dloc, d_global, w, d_collapse, d_grav, d_ref, floor):
    """Gamma = (super-critical fraction) x (does the eraser propagate?).

    U: fraction of the ensemble whose LOCAL dimension exceeds the bound-state
       threshold, softened by w (the model's one soft parameter). Because the
       ensemble has dispersion, U > 0 for any mean below d_collapse: the decay
       does NOT switch itself off at 4. This is why the flow does not stall at
       the upper threshold -- and it is a genuine prediction of the dispersion,
       not something put in by hand.
    G: the propagating-degree-of-freedom count of the long-range interaction.
       'dof'    -> D(D-3)/2, the linearised-GR graviton count in D spacetime
                   dimensions. Its zero at D = 3 is a theorem, not a dial.
       'linear' -> (d - d_grav), a free threshold. Used to show what happens
                   to the halt when the zero is moved.
    """
    u = float(np.mean(1.0 / (1.0 + np.exp(-(dloc - d_collapse) / w))))
    if floor == 'dof':
        g = max(0.0, d_global * (d_global - 3.0) / 2.0)
        g_ref = max(1e-12, d_ref * (d_ref - 3.0) / 2.0)
    else:
        g = max(0.0, d_global - d_grav)
        g_ref = max(1e-12, d_ref - d_grav)
    return u * g / g_ref, u, g / g_ref


# --------------------------------------------------------------------------
# the move
# --------------------------------------------------------------------------


def contract(x, eta, theta=THETA):
    """Erase directional support: shrink the thinnest still-resolvable axis.

    "Still resolvable" = principal std above theta x the IR extent, i.e. the
    direction is visible inside the measurement ball. Directions die
    thinnest-first, so the measured dimension falls smoothly instead of in unit
    jumps. The move knows nothing about any target: applied often enough it
    takes the cloud to a line and the estimator to ~1.
    """
    mu = x.mean(0)
    c = x - mu
    val, vec = np.linalg.eigh(c.T @ c / len(c))
    std = np.sqrt(np.maximum(val, 0.0))
    active = np.where(std > theta * std.max())[0]
    if active.size == 0:
        return x, None
    v = vec[:, active[0]]
    return (c - np.outer((c @ v) * eta, v)) + mu, float(std[active[0]])


# --------------------------------------------------------------------------
# the flow
# --------------------------------------------------------------------------


def asymptote(traj, max_tail=120):
    """Fit d(step) = d_inf + A * exp(-step/tau) over the tail of the flow.

    A cross-check on the halting value, not the primary number. The primary
    number is the dimension at which the rate actually fell below the stopping
    tolerance; this fit exists to show that the flow was still converging to
    the same place rather than crawling past it.

    Two guards, both learned the hard way on the first sweep: the tail adapts
    to the trajectory length (a fixed 120-step tail silently returned nothing
    for every run shorter than 140 steps, which quietly deleted the fastest
    half of the floor block from its own regression), and tau is capped at the
    tail length (an unbounded grid fits a huge slow mode to a nearly flat tail
    and extrapolates to a fixed point far below anything the run visited --
    that is what produced d_inf = 2.305 for a run that ended at 3.111).
    """
    if len(traj) < 50:
        return None
    tail = int(min(max_tail, max(40, len(traj) // 3)))
    s = np.array([p['step'] for p in traj[-tail:]], float)
    d = np.array([p['d'] for p in traj[-tail:]], float)
    best = None
    for tau in np.geomspace(3.0, float(tail), 80):
        x = np.exp(-(s - s[0]) / tau)
        a = np.vstack([np.ones_like(x), x]).T
        coef, *_ = np.linalg.lstsq(a, d, rcond=None)
        rms = float(np.sqrt(np.mean((a @ coef - d) ** 2)))
        if best is None or rms < best['rms']:
            best = dict(d_inf=float(coef[0]), amp=float(coef[1]),
                        tau=float(tau), rms=rms, tail=tail)
    return best


def run_flow(cal, seed=0, n=3000, ambient=AMBIENT, decay=0.794, w=0.50,
             eta0=0.60, d_collapse=4.0, d_grav=3.0, floor='dof',
             gamma_stop=1e-3, max_steps=400, nprobe=400, quiet=True,
             f_r=F_R, no_dispersion=False, capture_at=()):
    """no_dispersion: replace every per-probe value by the ensemble mean, so
    the rate sees a delta-function ensemble. Isolates how much of the flow's
    ability to pass through the upper threshold comes from real spread and how
    much from the softening width w alone.
    capture_at: dimensions at which to keep a copy of the state, for later
    scale-resolved inspection."""
    rng = np.random.default_rng(seed)
    x = init_state(rng, n, ambient, decay)
    d0 = apply_cal(measure(x, rng, nprobe, f_r)[0], cal)
    traj, halt, captures, pending = [], 'max_steps', {}, sorted(capture_at)[::-1]
    for step in range(max_steps):
        raw, dloc_raw, _ = measure(x, rng, nprobe, f_r)
        if not np.isfinite(raw):
            halt = 'estimator_starved'
            break
        d = apply_cal(raw, cal)
        dloc = apply_cal(dloc_raw, cal)
        if no_dispersion:
            dloc = np.full_like(dloc, float(np.mean(dloc)))
        while pending and d <= pending[0]:
            captures[pending.pop(0)] = x.copy()
        gam, u, g = instability_rate(dloc, d, w, d_collapse, d_grav, d0, floor)
        traj.append(dict(step=step, d=float(d), gamma=float(gam), u=float(u),
                         g=float(g), disp=float(np.std(dloc))))
        if not quiet and step % 40 == 0:
            print(f"  step {step:4d}  d = {d:6.3f}  Gamma = {gam:8.5f}  "
                  f"(U = {u:5.3f}, G = {g:5.3f})")
        if gam < gamma_stop:
            halt = 'rate_vanished'
            break
        x, _ = contract(x, eta0 * gam)
    d_final = apply_cal(measure(x, rng, nprobe, f_r)[0], cal)
    return dict(params=dict(seed=seed, n=n, ambient=ambient, decay=decay, w=w,
                            eta0=eta0, d_collapse=d_collapse, d_grav=d_grav,
                            floor=floor, gamma_stop=gamma_stop,
                            max_steps=max_steps, nprobe=nprobe, f_r=f_r,
                            theta=THETA, no_dispersion=no_dispersion),
                captures=captures,
                d_initial=float(d0), d_final=float(d_final), steps=len(traj),
                halt=halt, fixed_point=asymptote(traj), trajectory=traj), x


# --------------------------------------------------------------------------
# the frozen state: dimension vs scale (the CDT-facing prediction)
# --------------------------------------------------------------------------


def resample_like(x, n, rng):
    """Draw n fresh events from the frozen state's covariance.

    Every move in the flow is a linear contraction of a Gaussian, so the frozen
    state is exactly Gaussian and this resampling is exact, not an
    approximation. It exists so the scale-resolved measurements below are not
    limited by the N the flow was run at: at N = 3000 the diffusion runs out of
    lattice before it runs out of scales, and any UV claim would be a
    finite-size statement.
    """
    c = x - x.mean(0)
    val, vec = np.linalg.eigh(c.T @ c / len(c))
    return (rng.standard_normal((n, x.shape[1])) *
            np.sqrt(np.maximum(val, 0.0))) @ vec.T


def isotropic_control(k, n, ambient, rng):
    """A cloud of exactly known dimension k, for the same instrument to read."""
    x = np.zeros((n, ambient))
    x[:, :k] = rng.standard_normal((n, k))
    return x


def dim_vs_scale(x, cal, rng, fracs=(0.35, 0.5, 0.7, 1.0, 1.4, 2.0),
                 nprobe=400):
    """Measured dimension as a function of the size of the measuring ball."""
    out = []
    for f in fracs:
        d, _, r = measure(x, rng, nprobe, f_r=f)
        if np.isfinite(d):
            out.append(dict(f_r=f, r=float(r), d=float(apply_cal(d, cal))))
    return out


def spectral_dim(x, k=12, nsrc=64, taus=(4, 8, 16, 32, 64, 128, 256), seed=0):
    """Return-probability spectral dimension: d_s(tau) = -2 dlogP/dlogtau.

    This is the SAME observable CDT reports, obtained the same way -- a
    diffusion on the geometry, not a box-counting proxy -- so the disagreement
    with CDT's short-distance value is a comparison of like with like. Small
    tau probes short distances (UV), large tau long distances (IR).
    """
    rng = np.random.default_rng(seed)
    n = len(x)
    _, idx = cKDTree(x).query(x, k=k + 1)
    rows = np.repeat(np.arange(n), k)
    a = csr_matrix((np.ones(len(rows)), (rows, idx[:, 1:].ravel())),
                   shape=(n, n))
    a = a.maximum(a.T)                                   # reversible walk
    deg = np.asarray(a.sum(1)).ravel()
    deg[deg == 0] = 1.0
    w = csr_matrix((1.0 / deg, (np.arange(n), np.arange(n)))) @ a
    src = rng.choice(n, size=min(nsrc, n), replace=False)
    p = np.zeros((n, len(src)))
    p[src, np.arange(len(src))] = 1.0
    ret, t = {}, 0
    for tau in sorted(taus):
        while t < tau:
            p = w.T @ p
            t += 1
        ret[tau] = float(np.mean(p[src, np.arange(len(src))]))
    ts = sorted(ret)
    return [dict(tau_lo=ts[i], tau_hi=ts[i + 1],
                 d_s=float(-2.0 * np.log(ret[ts[i + 1]] / ret[ts[i]]) /
                           np.log(ts[i + 1] / ts[i])))
            for i in range(len(ts) - 1)
            if ret[ts[i]] > 0 and ret[ts[i + 1]] > 0]


# --------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--n', type=int, default=3000)
    ap.add_argument('--w', type=float, default=0.50)
    ap.add_argument('--eta0', type=float, default=0.60)
    ap.add_argument('--decay', type=float, default=0.794)
    ap.add_argument('--d-collapse', type=float, default=4.0)
    ap.add_argument('--d-grav', type=float, default=3.0)
    ap.add_argument('--floor', default='dof', choices=('dof', 'linear'))
    ap.add_argument('--max-steps', type=int, default=1200)
    ap.add_argument('--n-frozen', type=int, default=30000)
    ap.add_argument('--out', default='dimension/results/flow_default.json')
    args = ap.parse_args()

    t0 = time.time()
    print("instrument calibration (measured dimension vs known dimension):")
    cal = calibration(args.n)
    for k, d in zip(cal['k'], cal['d_raw']):
        print(f"   true {k:.0f}  ->  raw {d:6.3f}")
    print(f"   fit: raw = {cal['slope']:.4f} * true + {cal['intercept']:+.4f}"
          f"   (max resid {cal['max_resid']:.3f})")

    print("\nflow:")
    r, x = run_flow(cal, seed=args.seed, n=args.n, w=args.w, eta0=args.eta0,
                    decay=args.decay, d_collapse=args.d_collapse,
                    d_grav=args.d_grav, floor=args.floor,
                    max_steps=args.max_steps, quiet=False)
    r['calibration'] = cal
    print(f"  d: {r['d_initial']:.3f} -> {r['d_final']:.3f} in {r['steps']} "
          f"steps ({r['halt']})")
    fp = r['fixed_point']
    if fp:
        print(f"  exponential-approach fit: d_inf = {fp['d_inf']:.3f}  "
              f"(tau = {fp['tau']:.0f} steps, rms {fp['rms']:.4f})")

    rng = np.random.default_rng(args.seed + 991)
    big = resample_like(x, args.n_frozen, rng)
    ctrl = isotropic_control(3, args.n_frozen, AMBIENT, rng)
    r['scale_profile'] = dim_vs_scale(big, cal, rng)
    r['scale_profile_control_d3'] = dim_vs_scale(ctrl, cal, rng)
    r['spectral_profile'] = spectral_dim(big, seed=args.seed)
    r['spectral_profile_control_d3'] = spectral_dim(ctrl, seed=args.seed)
    r['n_frozen'] = args.n_frozen

    print(f"\nfrozen state at N = {args.n_frozen} (control: exact d = 3 cloud, "
          "same instrument)")
    print("dimension vs ball radius (UV first):")
    for b, c in zip(r['scale_profile'], r['scale_profile_control_d3']):
        print(f"  R = {b['f_r']:.2f} x IR extent   d = {b['d']:6.3f}"
              f"   [control {c['d']:.3f}]")
    print("spectral dimension by diffusion time (UV first):")
    for b, c in zip(r['spectral_profile'], r['spectral_profile_control_d3']):
        print(f"  tau {b['tau_lo']:4d}-{b['tau_hi']:4d}   d_s = {b['d_s']:6.3f}"
              f"   [control {c['d_s']:.3f}]")

    r['wall_seconds'] = round(time.time() - t0, 1)
    with open(args.out, 'w') as f:
        json.dump(r, f, indent=1)
    print(f"\nwrote {args.out}  ({r['wall_seconds']} s)")


if __name__ == '__main__':
    main()

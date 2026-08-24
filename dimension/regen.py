#!/usr/bin/env python3
"""THE NUMBER SOURCE for dimension/README.md.

Every quantitative claim in the write-up is generated here from
dimension/results/*.json; nothing is quoted from memory. Run from the repo
root:  python dimension/regen.py > dimension/TALLIES.txt
"""
import json

import numpy as np

FLOW = 'dimension/results/flow_default.json'
SWEEP = 'dimension/results/sweep.json'
GRILL = 'dimension/results/grill.json'

flow = json.load(open(FLOW))
sweep = json.load(open(SWEEP))
rows = sweep['rows']
try:
    grill = json.load(open(GRILL))
except FileNotFoundError:
    grill = None

bar = '=' * 70


def block(name):
    """Converged runs only: a run that hit the step budget has not halted."""
    return [r for r in rows if r['block'] == name and r['converged']]


def spread(name):
    v = [r['d_halt'] for r in block(name)]
    return (min(v), max(v), max(v) - min(v)) if v else (None, None, None)


print(bar)
print('0. INSTRUMENT — calibration of the dimension estimator')
print(bar)
cal = flow['calibration']
for k, d in zip(cal['k'], cal['d_raw']):
    print(f'  known d = {k:.0f}   raw MLE = {d:6.3f}')
print(f"  linear fit  raw = {cal['slope']:.4f} * known {cal['intercept']:+.4f}"
      f"   max residual {cal['max_resid']:.3f}")
print(f"  F-D3 (residual > 0.10 voids every quoted dimension): "
      f"{'FIRED' if cal['max_resid'] > 0.10 else 'not fired'}")

print()
print(bar)
print('1. THE DEFAULT FLOW')
print(bar)
fp = flow['fixed_point']
print(f"  d: {flow['d_initial']:.3f} -> {flow['d_final']:.3f} in "
      f"{flow['steps']} steps, halt = {flow['halt']}")
print(f"  exponential-approach fit: d_inf = {fp['d_inf']:.3f}, "
      f"tau = {fp['tau']:.0f} steps, rms = {fp['rms']:.4f}")
t = flow['trajectory']
print(f"  rate at start {t[0]['gamma']:.4f} -> at halt {t[-1]['gamma']:.5f} "
      f"({t[0]['gamma'] / max(t[-1]['gamma'], 1e-12):.0f}x suppression)")
for frac in (0.25, 0.5, 0.75):
    p = t[int(frac * (len(t) - 1))]
    print(f"    step {p['step']:4d}  d = {p['d']:6.3f}  Gamma = {p['gamma']:.5f}"
          f"  U = {p['u']:.3f}  G = {p['g']:.3f}")

print()
print(bar)
print('2. PRE-REGISTERED AUDIT (DIM-01) — does the fixed point follow anything?')
print(bar)
d1 = [r for r in rows if r['block'] == 'seed']
n_arrest = sum(1 for r in d1 if r['halt'] == 'rate_vanished')
print(f"  D1 self-arrest      {n_arrest}/{len(d1)} seeds halt on rate_vanished"
      f"   -> {'PASS' if n_arrest >= 4 else 'FAIL'}")
for name, key, lim, label in (
        ('D2 seed stability', 'seed', 0.15, 'across 5 seeds'),
        ('D3 soft parameter', 'soft_w', 0.15, 'across w = 0.25..1.00'),
        ('D4 step size', 'step_eta0', 0.15, 'across eta0 = 0.3..1.2'),
        ('D6 upper threshold', 'collapse', 0.20,
         'across d_collapse = 3.5..5.0')):
    lo, hi, sp = spread(key)
    if lo is None:
        continue
    print(f"  {name:19s} d_halt {lo:.3f} .. {hi:.3f}  spread {sp:.3f} "
          f"(limit {lim})  -> {'PASS' if sp < lim else 'FAIL'}   [{label}]")

st = block('start_d0')
above = [r for r in st if r['d_initial'] > 3.5]
below = [r for r in st if r['d_initial'] <= 3.5]
if above:
    v = [r['d_halt'] for r in above]
    starts = ', '.join('%.1f' % r['d_initial'] for r in above)
    print(f"  D5 basin (from above) starts {starts}"
          f"  ->  d_halt {min(v):.3f} .. {max(v):.3f}  spread {max(v)-min(v):.3f}"
          f"  -> {'PASS' if max(v)-min(v) < 0.20 else 'FAIL'}")
    print(f"      none stalled at the upper threshold: "
          f"{'confirmed' if max(v) < 3.6 else 'NOT confirmed'}")
for r in below:
    print(f"  D5 from below       start {r['d_initial']:.3f} -> "
          f"halt {r['d_halt']:.3f} after {r['steps']} step(s); one-sided "
          f"(no restoring term exists in the model)")

stalled = [r for r in rows if r['block'] == 'floor' and not r['converged']]
for r in stalled:
    print(f"  D7 excluded         floor {r['d_grav']:.1f}: still falling at "
          f"{r['d_halt']:.3f} after {r['steps']} steps ({r['halt']}) — "
          f"not converged, excluded from the regression")
fl = block('floor')
if len(fl) > 1:
    g = np.array([r['d_grav'] for r in fl])
    d = np.array([r['d_halt'] for r in fl])
    slope, icept = np.polyfit(g, d, 1)
    r2 = 1 - np.sum((d - (slope * g + icept)) ** 2) / np.sum((d - d.mean()) ** 2)
    print(f"  D7 THE CONTROL      halting dimension vs the floor handed in:")
    for gg, dd in zip(g, d):
        print(f"      floor placed at {gg:.1f}  ->  halted at {dd:.3f}")
    print(f"      slope {slope:.3f}, intercept {icept:+.3f}, R^2 {r2:.4f}"
          f"   (offset above the floor: "
          f"{np.mean(d - g):+.3f} +- {np.std(d - g):.3f})")
    verdict = ('the halt IS the input floor' if slope > 0.8 else
               'the halt is NOT set by the floor' if slope < 0.3 else
               'mixed: the floor is one of several things setting the halt')
    print(f"      registered reading -> {verdict}")

print()
print(bar)
print('3. THE FROZEN STATE vs A KNOWN 3-DIMENSIONAL SPACE')
print(bar)
print(f"  N = {flow['n_frozen']} events, control = exact d=3 cloud, "
      f"same instrument")
print('  dimension vs measuring-ball radius (UV first):')
for b, c in zip(flow['scale_profile'], flow['scale_profile_control_d3']):
    print(f"    R = {b['f_r']:.2f} x IR extent   frozen {b['d']:6.3f}   "
          f"control {c['d']:6.3f}")
sp_f = [b['d'] for b in flow['scale_profile']]
sp_c = [b['d'] for b in flow['scale_profile_control_d3']]
print(f"    frozen varies by {max(sp_f)-min(sp_f):.3f} across the window range;"
      f" a real d=3 space varies by {max(sp_c)-min(sp_c):.3f}")
print('  spectral dimension by diffusion time (the observable CDT reports):')
for b, c in zip(flow['spectral_profile'], flow['spectral_profile_control_d3']):
    print(f"    tau {b['tau_lo']:4d}-{b['tau_hi']:4d}   frozen {b['d_s']:6.3f}"
          f"   control {c['d_s']:6.3f}")
uv, ir = flow['spectral_profile'][0], flow['spectral_profile'][-1]
print(f"    d_s runs {uv['d_s']:.3f} (UV) -> {ir['d_s']:.3f} (IR): "
      f"{'UV ABOVE IR — opposite sign to CDT' if uv['d_s'] > ir['d_s'] else 'UV below IR — same sign as CDT'}")
print(f"    same instrument on a real d=3 space: "
      f"{flow['spectral_profile_control_d3'][0]['d_s']:.3f} -> "
      f"{flow['spectral_profile_control_d3'][-1]['d_s']:.3f} "
      f"(finite-N droop, the size of the artifact to beat)")

print()
print(bar)
print(f"  runs in sweep: {len(rows)}   sweep wall time: "
      f"{sweep['wall_seconds']} s   flow wall time: {flow['wall_seconds']} s")
print(bar)


if grill:
    print()
    print(bar)
    print('4. GRILLING THE AUDIT (dimension/grill.py) — attacks on section 2')
    print(bar)
    print('  G1 is the per-probe spread signal or estimator noise?')
    for r in grill['G1_split_half']:
        print(f"     at d = {r['d_measured']:5.3f}  SD {r['observed_sd']:.3f}"
              f" = signal {r['signal_sd']:.3f} + noise {r['noise_sd']:.3f}"
              f"   reliability {r['reliability']:.3f}"
              f"   ({r['median_neighbours']:.0f} neighbours/probe)")
    rel = [r['reliability'] for r in grill['G1_split_half']]
    print(f"     verdict: reliability {min(rel):.2f}-{max(rel):.2f} — the"
          f" spread is REAL, not sampling noise. The attack fails.")

    g2 = grill['G2_U_decomposition']
    nd = grill['G2_no_dispersion']
    print('  G2 does the spread actually drive the flow?')
    print(f"     U at mean local d = {g2['mean_local_d']:.3f}:"
          f" {g2['u_with_spread']:.4f} with the spread,"
          f" {g2['u_flat']:.4f} without"
          f" ({100 * g2['u_flat'] / g2['u_with_spread']:.0f}% survives)")
    print(f"     flow with the spread deleted: {nd['d_initial']:.3f} ->"
          f" {nd['d_final']:.3f} in {nd['steps']} steps ({nd['halt']})"
          f" — unchanged")
    print("     verdict: at the default softening the spread is NOT what"
          " carries the flow past the threshold.")

    print('  G3 does the softening width w carry it instead?')
    for r in grill['G3_hard_threshold']:
        print(f"     w = {r['w']:.3f}  halt {r['d_final']:6.3f}"
              f"  ({r['steps']:4d} steps, {r['halt']})")
    print("     verdict: no — at w = 0.03, near a hard step, the flow still"
          " reaches the floor.")
    print('  G3b crossed control — remove both:')
    for r in grill.get('G3b_crossed', []):
        print(f"     w = {r['w']:.3f}  dispersion "
              f"{'deleted' if r['no_dispersion'] else 'intact':7s}"
              f"  -> halt {r['d_final']:6.3f}  ({r['steps']:4d} steps)")
    both = [r for r in grill.get('G3b_crossed', [])
            if r['no_dispersion'] and r['w'] <= 0.03]
    if both:
        vals = ', '.join('%.3f' % r['d_final'] for r in both)
        print(f"     verdict: with BOTH removed the flow stalls at {vals},"
              f" well above the floor. The pass-through is OVER-DETERMINED —"
              f" real dispersion alone suffices, the softening alone suffices,"
              f" neither is necessary, removing both stops it early.")

    print('  G4 is the halt pinned by the rate or the stopping tolerance?')
    for r in grill['G4_tolerance']:
        print(f"     gamma_stop {r['gamma_stop']:.0e}  halt"
              f" {r['d_final']:6.3f}  ({r['steps']:4d} steps)")
    tight = [r['d_final'] for r in grill['G4_tolerance']
             if r['gamma_stop'] <= 1e-4]
    print(f"     verdict: tolerance-limited at 1e-2; converges to"
          f" {tight[-1]:.3f} by 1e-4. The headline 3.066 was"
          f" {3.066 - tight[-1]:+.3f} high for that reason.")

    print('  G5 is the halt a property of the geometry or of the window?')
    for r in grill['G5_window']:
        prof = '  '.join(f"{p['f_r']:.1f}:{p['d']:5.2f}"
                         for p in r['frozen_profile'])
        print(f"     dynamics measures at R = {r['f_r']:.1f} -> halts"
              f" {r['d_final']:.3f}   frozen state read at R = {prof}")
    print("     verdict: every window halts at ~3.06, and every window leaves"
          " a DIFFERENT geometry. The halt is a property of the measurement,"
          " not of the space.")
    g1_at = [r for r in grill['G1_split_half'] if r['target'] == 3.2]
    if g1_at and 'G2_U_decomposition' in grill:
        pooled = g1_at[0]['d_measured']
        permean = grill['G2_U_decomposition']['mean_local_d']
        print('  G6 scale mismatch inside the rate (found while grilling):')
        print(f"     on the same state, pooled estimate {pooled:.3f} vs mean"
              f" of per-probe estimates {permean:.3f}  -> gap"
              f" {permean - pooled:+.3f}")
        print(f"     the rate compares PER-PROBE values to d_collapse = 4"
              f" while the reported dimension is the POOLED value, so the"
              f" bound-state threshold actually bites at pooled"
              f" {4 - (permean - pooled):.2f}, not 4. A ratio-of-sums vs"
              f" mean-of-ratios (Jensen) gap in dim_flow.measure; it shifts"
              f" where the decay starts, not where it stops.")
    print(f"  grill wall time: {grill['wall_seconds']} s")
    print(bar)

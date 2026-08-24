#!/usr/bin/env python3
"""Robustness audit of the dimensional fixed point.

A single run that lands on 3 proves nothing on its own: the question is whether
3 is a property of the dynamics or a restatement of an input. Each block below
varies one thing and asks whether the halting dimension follows it.

The decisive block is FLOOR. If the fixed point tracks d_grav, then "the loop
chose three" means "the loop found the zero we handed it", and the emergent
content is the approach, not the number.

Run:  python dimension/dim_sweep.py --out dimension/results/sweep.json
"""
import argparse
import json
import time

import numpy as np

from dim_flow import calibration, run_flow

DEFAULTS = dict(n=3000, decay=0.794, w=0.50, eta0=0.60, d_collapse=4.0,
                d_grav=3.0, floor='dof', max_steps=900, nprobe=400)

# decay values calibrated in dim_flow to give these initial dimensions at n=3000
DECAY_FOR_D0 = {4.6: 0.740, 6.1: 0.794, 8.3: 0.850, 10.0: 0.890, 2.6: 0.580}

BLOCKS = [
    ('seed',        [dict(seed=s) for s in range(5)]),
    ('soft_w',      [dict(w=v, seed=1) for v in (0.25, 0.50, 0.75, 1.00)]),
    ('step_eta0',   [dict(eta0=v, seed=1) for v in (0.30, 0.60, 1.20)]),
    ('start_d0',    [dict(decay=DECAY_FOR_D0[v], seed=1)
                     for v in (2.6, 4.6, 6.1, 8.3, 10.0)]),
    ('size_n',      [dict(n=v, seed=1) for v in (2000, 3000, 5000)]),
    ('collapse',    [dict(d_collapse=v, seed=1, floor='linear', d_grav=3.0)
                     for v in (3.5, 4.0, 4.5, 5.0)]),
    ('floor',       [dict(floor='linear', d_grav=v, seed=1, max_steps=2500)
                     for v in (1.5, 2.0, 2.5, 3.0, 3.5)]),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='dimension/results/sweep.json')
    args = ap.parse_args()

    t0 = time.time()
    cals = {}
    rows = []
    for block, variants in BLOCKS:
        print(f"\n=== {block} " + "=" * (58 - len(block)))
        for v in variants:
            p = dict(DEFAULTS)
            p.update(v)
            n = p['n']
            if n not in cals:
                cals[n] = calibration(n)
            r, _ = run_flow(cals[n], **p)
            fp = r['fixed_point']
            row = dict(block=block, **{k: p[k] for k in
                                       ('seed', 'n', 'decay', 'w', 'eta0',
                                        'd_collapse', 'd_grav', 'floor')},
                       d_initial=r['d_initial'], d_halt=r['d_final'],
                       converged=(r['halt'] == 'rate_vanished'),
                       steps=r['steps'], halt=r['halt'],
                       d_inf=(fp or {}).get('d_inf'),
                       fit_rms=(fp or {}).get('rms'))
            rows.append(row)
            di = row['d_inf']
            print(f"  {block:9s} "
                  f"seed={p['seed']} n={n} w={p['w']:.2f} eta0={p['eta0']:.2f} "
                  f"dc={p['d_collapse']:.1f} dg={p['d_grav']:.1f} "
                  f"floor={p['floor']:6s} | "
                  f"d {r['d_initial']:6.3f} -> {r['d_final']:6.3f} "
                  f"in {r['steps']:4d} ({r['halt']}) "
                  f"d_inf={'  n/a' if di is None else f'{di:6.3f}'}")

    out = dict(rows=rows, defaults=DEFAULTS,
               calibrations={str(k): v for k, v in cals.items()},
               wall_seconds=round(time.time() - t0, 1))
    with open(args.out, 'w') as f:
        json.dump(out, f, indent=1)

    print("\n" + "=" * 70)
    print("SUMMARY — does the fixed point follow the varied quantity?")
    print("=" * 70)
    for block, _ in BLOCKS:
        vals = [r['d_halt'] for r in rows
                if r['block'] == block and r['converged']]
        if not vals:
            continue
        print(f"  {block:9s}  d_halt range {min(vals):6.3f} .. {max(vals):6.3f}"
              f"   spread {max(vals) - min(vals):.3f}   "
              f"({len(vals)} converged runs)")
    fl = [r for r in rows if r['block'] == 'floor' and r['converged']]
    if len(fl) > 1:
        g = np.array([r['d_grav'] for r in fl])
        d = np.array([r['d_halt'] for r in fl])
        slope, icept = np.polyfit(g, d, 1)
        print(f"\n  floor block: d_halt = {slope:.3f} * d_grav {icept:+.3f}"
              f"   (slope 1 => the halt IS the input floor)")
    print(f"\nwrote {args.out}  ({out['wall_seconds']} s)")


if __name__ == '__main__':
    main()

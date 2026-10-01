#!/usr/bin/env python3
"""ELIM-01 production: generate ensembles and store raw per-config measurements.

One process per ensemble; raw arrays go to elimination/raw/<tag>.npz and
nothing is analysed here, so the analysis (analyse.py) can be re-run without
regenerating configurations.

Run from repo root:  python elimination/run_elim.py
"""
import os
import sys
import time
from multiprocessing import Pool

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'contrast'))
sys.path.insert(0, HERE)
from su2_lattice import Lattice  # noqa: E402
from observables import (ape_smear, glueball_slices, polyakov,  # noqa: E402
                         wilson_rt)

SMEAR_LEVELS = (2, 6, 12, 20)        # cumulative APE steps, eps = 0.5
WILSON_LEVEL = 12
TMAX = 4
ENSEMBLES = [
    dict(tag='b2.3_L12', beta=2.3, L=12, meas=900, wilson_every=2),
    dict(tag='b2.4_L12', beta=2.4, L=12, meas=900, wilson_every=2),
    dict(tag='b2.3_L10', beta=2.3, L=10, meas=1300, wilson_every=1),
    dict(tag='b2.2_L8', beta=2.2, L=8, meas=2000, wilson_every=1),
    dict(tag='b2.3_L8', beta=2.3, L=8, meas=2000, wilson_every=1),
    dict(tag='b2.3_L6', beta=2.3, L=6, meas=3000, wilson_every=1),
]
THERM, SKIP = 300, 2


def run(e):
    t0 = time.time()
    lat = Lattice(4, e['L'], e['beta'], seed=int(1000 * e['beta']) + e['L'])
    for _ in range(THERM):
        lat.sweep()
    G, P, W, plaq = [], [], [], []
    for k in range(e['meas']):
        for _ in range(SKIP):
            lat.sweep()
        U = lat.U
        plaq.append(lat.plaquette())
        P.append(polyakov(U))
        S, done, ops, Sw = U, 0, [], None
        for n in SMEAR_LEVELS:
            S = ape_smear(S, n - done)
            done = n
            ops.append(glueball_slices(S))
            if n == WILSON_LEVEL:
                Sw = S
        G.append(ops)
        if k % e['wilson_every'] == 0:
            W.append(wilson_rt(U, Sw, e['L'] // 2, TMAX))
        if k % 100 == 0:
            print(f"  {e['tag']}: {k}/{e['meas']}  "
                  f"({time.time() - t0:.0f}s)", flush=True)
    np.savez_compressed(os.path.join(HERE, 'raw', e['tag'] + '.npz'),
                        G=np.array(G), P=np.array(P), W=np.array(W),
                        plaq=np.array(plaq), beta=e['beta'], L=e['L'],
                        smear_levels=np.array(SMEAR_LEVELS),
                        wilson_level=WILSON_LEVEL, therm=THERM, skip=SKIP,
                        wilson_every=e['wilson_every'],
                        seconds=time.time() - t0)
    print(f"DONE {e['tag']} in {time.time() - t0:.0f}s", flush=True)
    return e['tag']


if __name__ == '__main__':
    os.makedirs(os.path.join(HERE, 'raw'), exist_ok=True)
    only = sys.argv[1:]
    todo = [e for e in ENSEMBLES if not only or e['tag'] in only]
    with Pool(4) as pool:
        pool.map(run, todo, chunksize=1)

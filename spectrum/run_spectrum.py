#!/usr/bin/env python3
"""SPEC-01 production: all-channel glueball correlators at beta 2.3, L 10.

Four independent streams run in parallel. Each stores per-bin sufficient
statistics per channel (sums of slice-averaged operators and of
<Y(t0) Y(t0+t)^T>), enough for every analysis SPEC.md allows, at a few MB
instead of tens.

Run from repo root:  python spectrum/run_spectrum.py
"""
import os
import sys
import time
from multiprocessing import Pool

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'contrast'))
sys.path.insert(0, os.path.join(HERE, '..', 'elimination'))
sys.path.insert(0, HERE)
from su2_lattice import Lattice, qdag, qmul  # noqa: E402
from observables import ape_smear  # noqa: E402
from cubic import CHANNELS, group, isotypic_basis, orbit  # noqa: E402

BETA, L = 2.3, 10
STREAMS, MEAS, THERM, SKIP, NBIN = 4, 1000, 300, 2, 5
LEVELS = (5, 15)
PROTOTYPES = [
    (-3, -1, 3, 1),                       # plaquette
    (-3, -3, -1, 3, 3, 1),                # 2x1 rectangle
    (-3, -2, -1, 2, 3, 1),                # 6-link bent ("chair")
    (-3, -2, -1, 3, 2, 1),                # 6-link twisted
    (-3, -3, -2, -1, 3, 1, 3, 2),         # 8-link chiral, orbit 48
    (-3, -3, -2, -1, 3, 2, 3, 1),         # 8-link chiral, orbit 48
]
T = L // 2 + 1


def build_projection():
    """Per channel: list over prototypes of (orbit elements, basis matrix)."""
    G = group()
    orbits = [orbit(p, G) for p in PROTOTYPES]
    proj = {}
    for ch in CHANNELS:
        proj[ch] = [(elems, isotypic_basis(perms, G, *ch))
                    for elems, perms in orbits]
    return orbits, proj


def loop_slices(S, path, cache):
    """sum over spatial x of (1/2)Tr(loop starting at x) per time slice."""
    d = np.zeros(3, int)
    P = None
    for s in path:
        i = abs(s)
        if s < 0:
            d[i - 1] -= 1
        key = (i, tuple(d))
        if key not in cache:
            cache[key] = np.roll(S[i], tuple(-d), axis=(1, 2, 3))
        link = cache[key] if s > 0 else qdag(cache[key])
        if s > 0:
            d[i - 1] += 1
        P = link if P is None else qmul(P, link)
    assert not d.any(), 'path is not closed'
    return P[..., 0].sum(axis=(1, 2, 3))


def stream(k):
    t0 = time.time()
    orbits, proj = build_projection()
    elems_all = [e for es, _ in orbits for e in es]
    lat = Lattice(4, L, BETA, seed=7000 + k)
    for _ in range(THERM):
        lat.sweep()
    acc = {ch: None for ch in CHANNELS}
    per_bin = MEAS // NBIN
    for m in range(MEAS):
        for _ in range(SKIP):
            lat.sweep()
        S, done, vals = lat.U, 0, []
        for n in LEVELS:
            S = ape_smear(S, n - done)
            done = n
            cache = {}
            vals.append({e: loop_slices(S, e, cache) for e in elems_all})
        b = m // per_bin
        for ch in CHANNELS:
            rows = []
            for lev in vals:
                for elems, B in proj[ch]:
                    if B.shape[1] == 0:
                        continue
                    O = np.array([lev[e] for e in elems])        # (n_p, L)
                    rows.append(B.T @ O)                          # (k, L)
            Y = np.concatenate(rows, axis=0)                      # (K, L)
            if acc[ch] is None:
                K = Y.shape[0]
                acc[ch] = dict(S1=np.zeros((NBIN, K)),
                               S2=np.zeros((NBIN, K, K, T)),
                               n=np.zeros(NBIN))
            a = acc[ch]
            a['S1'][b] += Y.mean(1)
            for t in range(T):
                a['S2'][b, :, :, t] += Y @ np.roll(Y, -t, axis=1).T / L
            a['n'][b] += 1
        if m % 100 == 0:
            print(f'  stream {k}: {m}/{MEAS}  ({time.time() - t0:.0f}s)',
                  flush=True)
    out = {}
    for (r, p), a in acc.items():
        tag = f"{r}{'p' if p > 0 else 'm'}"
        for key, v in a.items():
            out[f'{tag}_{key}'] = v
    np.savez_compressed(os.path.join(HERE, 'raw', f'stream{k}.npz'),
                        beta=BETA, L=L, levels=np.array(LEVELS), skip=SKIP,
                        therm=THERM, seconds=time.time() - t0, **out)
    print(f'DONE stream {k} in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    os.makedirs(os.path.join(HERE, 'raw'), exist_ok=True)
    with Pool(STREAMS) as pool:
        pool.map(stream, range(STREAMS), chunksize=1)

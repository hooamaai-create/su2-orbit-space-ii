#!/usr/bin/env python3
"""MIRROR-01 production: 15 ensembles, per-bin sufficient statistics.

Run from repo root:  python mirror/run_mirror.py
"""
import os
import sys
import time
from multiprocessing import Pool

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'contrast'))
from open_lattice import OpenLattice, ape_smear_open  # noqa: E402
from su2_lattice import qmul  # noqa: E402

BETA, LT, LXS, DS = 2.3, 16, (6, 8, 12), (1, 2, 3, 4, 6)
MEAS, THERM, SKIP, NBIN, LEVELS = 3000, 300, 2, 20, (3, 10)
T = LT // 2 + 1


def slice_plaq(lat, S, mu, nu):
    v = lat.pvalid[(mu, nu)]
    if not v.any():
        return None
    P = qmul(qmul(S[mu], lat.sh(S[nu], mu, 1)),
             qmul(lat.qdag_(lat.sh(S[mu], nu, 1)), lat.qdag_(S[nu])))[..., 0]
    return (P * v).reshape(LT, -1).sum(1)


def operators(lat):
    U = lat.U
    ops = {}
    if lat.shape[2] == 1:
        ops['tx'] = slice_plaq(lat, U, 0, 1)
    else:
        S, done = U, 0
        for n in LEVELS:
            S = ape_smear_open(lat, S, n - done)
            done = n
            mix = slice_plaq(lat, S, 1, 2) + slice_plaq(lat, S, 1, 3)
            ops[f'mix{n}'] = mix
            ops[f'trv{n}'] = slice_plaq(lat, S, 2, 3)
    P = U[1][:, 0].copy()
    for x in range(1, lat.shape[1]):
        P = qmul(P, U[1][:, x])
    ops['poly'] = P[..., 0].reshape(LT, -1).mean(1)
    return ops


def run(job):
    Lx, d = job
    t0 = time.time()
    lat = OpenLattice(BETA, (LT, Lx, d, d), (2, 3), seed=31 * Lx + d)
    from su2_lattice import qdag
    lat.qdag_ = qdag
    for _ in range(THERM):
        lat.sweep()
    names, S1, S2, cnt = None, None, None, np.zeros(NBIN)
    per = MEAS // NBIN
    for m in range(MEAS):
        for _ in range(SKIP):
            lat.sweep()
        ops = operators(lat)
        if names is None:
            names = sorted(ops)
            K = len(names)
            S1 = np.zeros((NBIN, K))
            S2 = np.zeros((NBIN, K, K, T))
        Y = np.array([ops[k] for k in names])
        b = m // per
        S1[b] += Y.mean(1)
        for t in range(T):
            S2[b, :, :, t] += Y @ np.roll(Y, -t, axis=1).T / LT
        cnt[b] += 1
    np.savez_compressed(os.path.join(HERE, 'raw', f'Lx{Lx}_d{d}.npz'),
                        names=np.array(names), S1=S1, S2=S2, n=cnt, Lx=Lx, d=d,
                        beta=BETA, LT=LT, seconds=time.time() - t0)
    print(f'DONE Lx={Lx} d={d} in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    os.makedirs(os.path.join(HERE, 'raw'), exist_ok=True)
    jobs = sorted([(Lx, d) for Lx in LXS for d in DS],
                  key=lambda j: -j[0] * j[1] ** 2)
    with Pool(4) as pool:
        pool.map(run, jobs, chunksize=1)

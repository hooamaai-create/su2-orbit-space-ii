#!/usr/bin/env python3
"""TWOPLANES-01 production. Run from repo root: python twoplanes/run_twoplanes.py"""
import os
import sys
import time
from multiprocessing import Pool

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
for p in ('', '../contrast', '../elimination'):
    sys.path.insert(0, os.path.join(HERE, p))
from plane_lattice import PlaneLattice  # noqa: E402
from observables import ape_smear  # noqa: E402
from su2_lattice import qdag, qmul  # noqa: E402

BETA, LT, LS, KAPPAS = 2.3, 16, (6, 8, 10), (0.0, 0.25, 0.5, 0.75, 1.0)
MEAS, THERM, SKIP, NBIN, LEVELS = 1500, 300, 2, 20, (3, 10)
T = LT // 2 + 1


def slice_plaq(S, mu, nu):
    P = qmul(qmul(S[mu], np.roll(S[nu], -1, axis=mu)),
             qmul(qdag(np.roll(S[mu], -1, axis=nu)), qdag(S[nu])))[..., 0]
    return P.reshape(LT, -1).sum(1)


def operators(U):
    ops, S, done = {}, U, 0
    for n in LEVELS:
        S = ape_smear(S, n - done)
        done = n
        ops[f'mix{n}'] = slice_plaq(S, 1, 2) + slice_plaq(S, 1, 3)
        ops[f'shB{n}'] = slice_plaq(S, 2, 3)
    P = U[1][:, 0].copy()
    for x in range(1, U.shape[2]):
        P = qmul(P, U[1][:, x])
    ops['poly'] = P[..., 0].reshape(LT, -1).mean(1)
    return ops


def run(job):
    L, kappa = job
    t0 = time.time()
    lat = PlaneLattice(BETA, (LT, L, L, L), kappa, seed=int(1000 * kappa) + L)
    for _ in range(THERM):
        lat.sweep()
    names = S1 = S2 = None
    n = np.zeros(NBIN)
    plaq = np.zeros((NBIN, 3))
    per = MEAS // NBIN
    for m in range(MEAS):
        for _ in range(SKIP):
            lat.sweep()
        ops = operators(lat.U)
        if names is None:
            names = sorted(ops)
            K = len(names)
            S1, S2 = np.zeros((NBIN, K)), np.zeros((NBIN, K, K, T))
        Y = np.array([ops[k] for k in names])
        b = m // per
        S1[b] += Y.mean(1)
        for t in range(T):
            S2[b, :, :, t] += Y @ np.roll(Y, -t, axis=1).T / LT
        n[b] += 1
        plaq[b] += [lat.plane_plaquette(0, 1), lat.plane_plaquette(2, 3),
                    lat.plane_plaquette(1, 2)]
    np.savez_compressed(os.path.join(HERE, 'raw', f'L{L}_k{kappa:.2f}.npz'),
                        names=np.array(names), S1=S1, S2=S2, n=n, plaq=plaq,
                        L=L, kappa=kappa, beta=BETA, LT=LT,
                        seconds=time.time() - t0)
    print(f'DONE L={L} kappa={kappa:.2f} in {time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    os.makedirs(os.path.join(HERE, 'raw'), exist_ok=True)
    jobs = sorted([(L, k) for L in LS for k in KAPPAS], key=lambda j: -j[0])
    with Pool(4) as pool:
        pool.map(run, jobs, chunksize=1)

#!/usr/bin/env python3
"""POST-HOC (not registered): the M2 particle-vs-flux-line fit using m(0->1).

Written after reading REPORT.txt, because at d = 2 and 3 the registered
statistic m(1->2) has errors of 0.3-0.9, too large to tell a constant from a
line, while m(0->1) is precise. m(0->1) is an upper bound on the lightest
energy in the channel (excited states raise it), but whether it stays fixed or
grows with the box is still the particle-vs-flux question. The d = 1 control is
included: there the channel is empty and m(0->1) is noise.
"""
import json
import os

import numpy as np
from scipy.stats import chi2 as chi2dist

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'results.json')))
LXS = (6, 8, 12)
print('POST-HOC: constant (particle) vs E = k*L (flux line), using m(0->1)')
out = {}
for d in (1, 2, 3, 4, 6):
    E = np.array([R[f'{L}_{d}']['glue']['m01'] for L in LXS])
    e = np.array([R[f'{L}_{d}']['glue']['m01_err'] for L in LXS])
    F = np.array([R[f'{L}_{d}']['flux']['m01'] for L in LXS])
    L = np.array(LXS, float)
    w = 1 / e ** 2
    c0 = np.sum(w * E) / np.sum(w)
    chi_c = float(np.sum(w * (E - c0) ** 2))
    k = np.sum(w * E * L) / np.sum(w * L ** 2)
    chi_f = float(np.sum(w * (E - k * L) ** 2))
    p = float(1 - chi2dist.cdf(chi_c, 2))
    tag = ('particle-like' if p > 0.01 and chi_f - chi_c > 9 else 'undecided')
    out[d] = dict(m01=E.tolist(), err=e.tolist(), const=float(c0), chi_const=chi_c,
                  p_const=p, chi_flux=chi_f, reading=tag)
    fl = ', '.join('n/a' if not np.isfinite(x) else f'{x:.2f}' for x in F)
    print(f'  d={d}: m(0->1) = ' + ', '.join(f'{a:.3f}({b * 1e3:.0f})' for a, b in zip(E, e))
          + f'   const {c0:.3f}, chi2 {chi_c:.1f}/2 (p={p:.2f})   flux-fit chi2 {chi_f:.1f}/2'
          f'   [flux line itself: {fl}]  -> {tag}')
json.dump(out, open(os.path.join(HERE, 'posthoc_m01.json'), 'w'), indent=1)

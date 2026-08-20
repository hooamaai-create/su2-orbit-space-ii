"""EXPLORATORY (run-14 data, post-hoc, labeled): per-config Gribov fits in the
dead ensemble D1 -> gamma_i, Z_i per configuration -> regress on lambda_i.
Also sliding-window fits (stabler than per-config) and candidate formulas."""
import json
import numpy as np

R = json.load(open('results14.json'))
e = R['ensembles']['D1']
lam = np.array(e['cfg']['lam'])
D = np.array(e['cfg']['D_shells'])          # (64, 4)
ph2 = np.array(e['phat2'])
n = len(lam)

# ---- per-config M2 fit: D_i(p) = Z * p2/(p2^2 + g^4); 4 points, 2 params ----
grid = np.linspace(1e-3, 1.5, 1500)
shapes = ph2[None, :] / (ph2[None, :] ** 2 + grid[:, None] ** 4)   # (G, 4)


def fit_one(d):
    # least squares over grid; equal weights (no per-config errors available)
    Z = (d[None, :] * shapes).sum(1) / (shapes ** 2).sum(1)
    c2 = ((d[None, :] - Z[:, None] * shapes) ** 2).sum(1)
    i = int(np.argmin(c2))
    # crude 1-sigma from delta-c2 = c2_min * (1 + 1/dof): use c2 doubling as band
    return grid[i], Z[i], c2[i]


gam = np.zeros(n); Z = np.zeros(n)
for i in range(n):
    gam[i], Z[i], _ = fit_one(D[i])

print(f"D1 dead ensemble, n={n}: per-config Gribov fits")
print(f"  gamma_i: mean {gam.mean():.3f}  std {gam.std(ddof=1):.3f}  "
      f"range [{gam.min():.3f}, {gam.max():.3f}]")
print(f"  Z_i:     mean {Z.mean():.2f}   std {Z.std(ddof=1):.2f}")


def rp(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a*b).sum() / np.sqrt((a*a).sum() * (b*b).sum()))


print(f"\n  corr(lambda, gamma_i) = {rp(lam, gam):+.3f}")
print(f"  corr(lambda, Z_i)     = {rp(lam, Z):+.3f}")
print(f"  corr(lambda, log gam) = {rp(lam, np.log(gam)):+.3f}")
print(f"  corr(sqrt(lam), gam)  = {rp(np.sqrt(lam), gam):+.3f}")
print(f"  corr(Z_i, gamma_i)    = {rp(Z, gam):+.3f}")

# ---- sliding windows of 16 configs sorted by lambda (stabler) ----
o = np.argsort(lam)
print("\n  sliding windows (16 configs, step 8), ensemble-style fit:")
print(f"  {'<lam>':8} {'gamma':7} {'Z':7}")
xs, gs, zs = [], [], []
for s in range(0, n - 15, 8):
    idx = o[s:s+16]
    Dm = D[idx].mean(0); De = D[idx].std(0, ddof=1) / 4.0
    w = 1 / De ** 2
    best = (1e9, 0, 0)
    for g in grid:
        sh = ph2 / (ph2 ** 2 + g ** 4)
        Zw = (w * Dm * sh).sum() / (w * sh ** 2).sum()
        c2 = ((Dm - Zw * sh) ** 2 * w).sum()
        if c2 < best[0]:
            best = (c2, g, Zw)
    xs.append(lam[idx].mean()); gs.append(best[1]); zs.append(best[2])
    print(f"  {xs[-1]:8.4f} {gs[-1]:7.3f} {zs[-1]:7.2f}")

xs, gs, zs = map(np.array, (xs, gs, zs))
# candidate formulas on windowed points
for name, X in (("gamma = a + b*lam", xs),
                ("gamma = a + b*sqrt(lam)", np.sqrt(xs)),
                ("gamma = a + b*log(lam)", np.log(xs))):
    A = np.vstack([np.ones_like(X), X]).T
    coef, res, *_ = np.linalg.lstsq(A, gs, rcond=None)
    pred = A @ coef
    ss = 1 - ((gs - pred) ** 2).sum() / ((gs - gs.mean()) ** 2).sum()
    print(f"  {name:26} a={coef[0]:+.3f} b={coef[1]:+.3f}  R2={ss:.3f}")
for name, X in (("Z = a + b*lam", xs), ("Z = a + b*log(lam)", np.log(xs))):
    A = np.vstack([np.ones_like(X), X]).T
    coef, *_ = np.linalg.lstsq(A, zs, rcond=None)
    pred = A @ coef
    ss = 1 - ((zs - pred) ** 2).sum() / ((zs - zs.mean()) ** 2).sum()
    print(f"  {name:26} a={coef[0]:+.3f} b={coef[1]:+.3f}  R2={ss:.3f}")

"""QUANT SWEEP (exploratory, labeled): mass & energy across runs 8-15.
Master table -> what drives IR suppression (mass)? what does energy do?"""
import json
import numpy as np

ASQS = {2.2: 0.561, 2.3: 0.4975, 2.4: 0.434, 2.5: 0.385, 2.6: 0.336}


def load(p):
    return json.load(open(p))


rows = []
for path, run in (('r12out/results12.json', 12), ('r13out/results13A.json', 13),
                  ('r13out/results13B.json', 13), ('r14out/results14.json', 14),
                  ('r15out/results15.json', 15)):
    R = load(path)
    for t, e in R['ensembles'].items():
        dp2 = [e['shells'][f'k2={s}']['D_mean'] * e['phat2'][s-1] for s in (1, 2, 3, 4)]
        de = [e['shells'][f'k2={s}']['D_err'] * e['phat2'][s-1] for s in (1, 2, 3, 4)]
        # mass proxy: s = ln(Dp2_shell2 / Dp2_shell1); >0 = IR-suppressed (massive)
        s = np.log(dp2[1] / dp2[0])
        s_err = np.sqrt((de[0]/dp2[0])**2 + (de[1]/dp2[1])**2)
        f = e['fits']
        d02 = f['M0']['chi2'] - f['M2']['chi2']
        d01 = f['M0']['chi2'] - f['M1']['chi2']
        gam = f['M2']['scale'] if d02 > 6 else np.nan
        m1 = f['M1']['scale'] if d01 > 6 else np.nan
        plaq = None
        if 'cfg' in e:
            plaq = float(np.mean(e['cfg']['plaq']))
        rows.append(dict(run=run, tag=t, beta=e['beta'], L=e['L'], n=e['n'],
                         ell=e['L']*ASQS[e['beta']], f1=e['f1_mean'],
                         lam=e['lam_mean'] if 'lam_mean' in e else float(np.mean(e['cfg']['lam'])),
                         s=s, s_err=s_err, gam=gam, m1=m1, Z=f['M0']['Z'],
                         plaq=plaq))

rows.sort(key=lambda r: r['f1'])
print(f"{'run':3} {'tag':4} {'b':4} {'L':3} {'ell':5} {'f1':6} {'lam':7} {'s(mass)':8} {'gam_lat':7} {'Z':6}")
for r in rows:
    gtxt = f"{r['gam']:.3f}" if np.isfinite(r['gam']) else "  --  "
    print(f"{r['run']:3} {r['tag']:4} {r['beta']:.1f} {r['L']:3} {r['ell']:5.2f} "
          f"{r['f1']:.3f} {r['lam']:7.4f} {r['s']:+8.3f} {gtxt:7} {r['Z']:6.2f}")

# ---- what drives the mass proxy s? univariate + rank correlations ----
def pear(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a - a.mean(); b = b - b.mean()
    return float((a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum()))


def rank(x):
    return np.argsort(np.argsort(x)).astype(float)


s = np.array([r['s'] for r in rows])
print(f"\nDRIVERS of mass proxy s (n={len(rows)} ensembles):  pearson / spearman")
for k in ('f1', 'ell', 'beta', 'L', 'lam', 'Z'):
    v = np.array([r[k] for r in rows], float)
    print(f"  s vs {k:4}: {pear(v, s):+.3f} / {pear(rank(v), rank(s)):+.3f}")

# two-variable horse race: f1 vs ell (partial correlations)
f1 = np.array([r['f1'] for r in rows]); ell = np.array([r['ell'] for r in rows])


def partial(x, y, z):
    rxy, rxz, ryz = pear(x, y), pear(x, z), pear(y, z)
    return (rxy - rxz*ryz)/np.sqrt((1-rxz**2)*(1-ryz**2))


print(f"\n  partial s~f1 | ell = {partial(f1, s, ell):+.3f}")
print(f"  partial s~ell | f1 = {partial(ell, s, f1):+.3f}")

# gamma structure where a scale exists
gs = [(r['beta'], r['L'], r['f1'], r['lam'], r['gam'], r['ell']) for r in rows if np.isfinite(r['gam'])]
print(f"\nGAMMA (lattice units) where M2 wins (n={len(gs)}):")
b = np.array([g[0] for g in gs]); Lg = np.array([g[1] for g in gs])
fg = np.array([g[2] for g in gs]); lg = np.array([g[3] for g in gs])
gg = np.array([g[4] for g in gs]); eg = np.array([g[5] for g in gs])
for k, v in (('beta', b), ('L', Lg), ('f1', fg), ('lam_mean', lg), ('ell', eg)):
    print(f"  gam vs {k:8}: pearson {pear(v, gg):+.3f}")
a = np.array([ASQS[x] for x in b])
print(f"  gam vs a(sqs) : pearson {pear(a, gg):+.3f}   <- lattice-spacing tie?")
# candidate: gam_lat = c * a^p  (log-log)
X = np.log(a); Y = np.log(gg)
p, c = np.polyfit(X, Y, 1)
pred = np.polyval([p, c], X)
r2 = 1 - ((Y-pred)**2).sum()/((Y-Y.mean())**2).sum()
print(f"  fit gam_lat = C * a^p:  p = {p:.2f}, C = {np.exp(c):.3f}, R2 = {r2:.3f}")
print(f"  -> gam_phys = gam_lat/a = C * a^(p-1); p=1 means gam_phys = const = {np.exp(c):.3f} sqrt(sigma)")

# ---- ENERGY: plaquette across the transition ----
print("\nENERGY (plaquette) across the mass transition:")
r8 = load('../su2-orbit-space/results/results8.json')['partO']
Ls = [k for k in r8 if k.isdigit()]
pl8 = [r8[k]['plaq'] for k in sorted(Ls, key=int)]
f8 = [r8[k]['f_ir'] for k in sorted(Ls, key=int)]
print(f"  beta=2.2 (run 8): L = {sorted(map(int, Ls))}")
print(f"    plaq = " + " ".join(f"{p:.4f}" for p in pl8))
print(f"    f_ir = " + " ".join(f"{f:.3f}" for f in f8))
print(f"    plaq rel. spread = {np.std(pl8)/np.mean(pl8)*100:.3f}%  while f_ir falls "
      f"{max(f8):.2f} -> {min(f8):.2f} and the mass turns ON")
pl_by_beta = {}
for r in rows:
    if r['plaq'] is not None:
        pl_by_beta.setdefault(r['beta'], []).append((r['L'], r['plaq'], r['s']))
for bb, lst in sorted(pl_by_beta.items()):
    for (L, p, sv) in sorted(lst):
        print(f"  beta={bb}, L={L}: plaq={p:.5f}  s={sv:+.3f}")

"""Run-15 verdicts: P6/P7 (H-II), P8 (curve), P9 (Z-law fork), P10 (localization),
P11 (energy floor), trans-horizon counts, and the M1-vs-M2 volume trend."""
import json
import numpy as np

R = json.load(open('results15.json'))


def pear(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum()))


def partial(x, y, z):
    rxy, rxz, ryz = pear(x, y), pear(x, z), pear(y, z)
    return (rxy - rxz*ryz)/np.sqrt((1-rxz**2)*(1-ryz**2))


def sig(r, n):
    return 0.5*np.log((1+r)/(1-r))*np.sqrt(max(n-4, 1))


def fit_M2_gamma(D, ph2, mask, w_from_err=True):
    Dm = D[mask].mean(0)
    De = D[mask].std(0, ddof=1)/np.sqrt(mask.sum())
    w = 1/De**2
    grid = np.linspace(1e-4, 2.0, 2000)
    best = (1e18, 0)
    cc = np.zeros(len(grid))
    for i, g in enumerate(grid):
        sh = ph2/(ph2**2 + g**4)
        Z = (w*Dm*sh).sum()/(w*sh**2).sum()
        cc[i] = ((Dm - Z*sh)**2*w).sum()
    i = int(cc.argmin()); ok = cc <= cc[i]+1
    return grid[i], grid[ok].min(), grid[ok].max()


print(f"gates pass: {R['gates']['all_pass']}\n")
for t, e in R['ensembles'].items():
    c = e['cfg']
    lam = np.array(c['lam']); f1 = np.array(c['f1'])
    pl = np.array(c['plaq']); D = np.array(c['D_shells'])
    ph2 = np.array(e['phat2']); n = e['n']
    Dk1 = D[:, 0]; Z = (D*ph2).mean(1)
    med = np.median(lam); near = lam <= med

    gN = fit_M2_gamma(D, ph2, near)
    gF = fit_M2_gamma(D, ph2, ~near)
    sep = "SEPARATED" if (gF[1] > gN[2] or gN[1] > gF[2]) else "overlap"
    p6 = "P6-PASS" if (gF[0] > gN[0] and gF[1] > gN[2]) else \
         ("ANTI" if (gN[0] > gF[0] and gN[1] > gF[2]) else "EMPTY/overlap")

    # P7 quartiles
    q = np.quantile(lam, [.25, .5, .75])
    masks = [lam <= q[0], (lam > q[0]) & (lam <= q[1]),
             (lam > q[1]) & (lam <= q[2]), lam > q[2]]
    gq = [fit_M2_gamma(D, ph2, m)[0] for m in masks]
    # Kendall tau of 4 points
    conc = sum(np.sign(gq[j]-gq[i]) for i in range(4) for j in range(i+1, 4))
    tau = conc/6.0

    rD = partial(lam, Dk1, pl); rZ = partial(lam, Z, pl)
    rf = partial(lam, f1, pl); rpl = pear(lam, pl)
    trans = int((lam < -1e-5).sum())

    print(f"=== {t} (beta={e['beta']}, L={e['L']}, n={n})  f1={e['f1_mean']:.3f}")
    print(f"  P6  gamma NEAR {gN[0]:.3f} [{gN[1]:.3f},{gN[2]:.3f}]  "
          f"FAR {gF[0]:.3f} [{gF[1]:.3f},{gF[2]:.3f}]  {sep}  -> {p6}")
    print(f"  P7  quartile gammas: " + " ".join(f"{g:.3f}" for g in gq) +
          f"   Kendall tau = {tau:+.2f}")
    print(f"  P8  r(lam,Dk1) = {rD:+.3f} ({sig(rD,n):+.1f}s)   [pred |r|<0.25]")
    print(f"  P9  r(lam,Z)   = {rZ:+.3f} ({sig(rZ,n):+.1f}s)   "
          f"[fork: (a) |r|<0.12  (b) in [-0.35,-0.15]]")
    print(f"  P10 r(lam,f1)  = {rf:+.3f} ({sig(rf,n):+.1f}s)   [pred > 0]")
    print(f"  P11 r(lam,plaq)= {rpl:+.3f} ({sig(rpl,n):+.1f}s)  [pred < 0]")
    print(f"  trans-horizon configs (lam < -1e-5): {trans}")
    f = e['fits']
    d01 = f['M0']['chi2']-f['M1']['chi2']; d02 = f['M0']['chi2']-f['M2']['chi2']
    d12 = f['M1']['chi2']-f['M2']['chi2']
    w = ('M1 (massive)' if d12 < -6 else ('M2 (Gribov)' if d12 > 6 else 'tie'))
    print(f"  ensemble fit: dchi2 M0-M1={d01:.1f} M0-M2={d02:.1f}  winner: {w}  "
          f"m={f['M1']['scale_phys']:.2f}sqs  gamma={f['M2']['scale_phys']:.2f}sqs\n")

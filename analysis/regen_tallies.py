"""THE NUMBER SOURCE. Every quantitative claim in README/description is generated
here from the raw results JSONs — nothing is quoted from memory. Run from repo
root: python analysis/regen_tallies.py > analysis/TALLIES.txt"""
import json
import glob
import numpy as np
from scipy.optimize import minimize

ASQS = {2.2: 0.561, 2.3: 0.4975, 2.4: 0.434, 2.5: 0.385, 2.6: 0.336}


def z(r):
    return np.arctanh(r)


def pear(a, b):
    a = np.asarray(a, float) - np.mean(a); b = np.asarray(b, float) - np.mean(b)
    return float((a*b).sum()/np.sqrt((a*a).sum()*(b*b).sum()))


def partial(x, y, zz):
    rxy, rxz, ryz = pear(x, y), pear(x, zz), pear(y, zz)
    return (rxy - rxz*ryz)/np.sqrt((1-rxz**2)*(1-ryz**2))


# ---------- 1. THE ACTIVATION CURVE: every blind (prediction, measurement) pair ----------
print("=" * 70)
print("1. ACTIVATION CURVE — all blind points (runs 11-19)")
print("=" * 70)
pairs = []
# run 11 (part I record, replicated here for the tally): predictions + partial r
for f in ('results/results11A.json', 'results/results11B.json'):
    try:
        R = json.load(open(f))
    except FileNotFoundError:
        continue
    for t, p in R['predictions'].items():
        e = R['ensembles'][t]
        pairs.append((f'r11:{t}', p['r_pred'], e['partial'], e['n']))
# runs 12-19 (this record): predictions + partial r at k2=1
for f in sorted(glob.glob('results/results1[2-9]*.json')):
    R = json.load(open(f))
    for t, p in R.get('predictions', {}).items():
        if t not in R['ensembles']:
            continue
        e = R['ensembles'][t]
        rm = e['shells']['k2=1']['partial'] if 'shells' in e else e['partial']
        pairs.append((f"{f.split('/')[-1][:-5]}:{t}", p['r_pred'], rm, e['n']))
chi2 = 0.0
worst = (0, '')
for tag, rp, rm, n in pairs:
    sz = 1/np.sqrt(n - 4)
    d = (z(rm) - z(rp))/sz
    chi2 += d*d
    if abs(d) > abs(worst[0]):
        worst = (d, tag)
    print(f"  {tag:22} pred {rp:+.3f}  meas {rm:+.3f}  dev {d:+5.2f} sigma  (n={n})")
print(f"\n  TOTAL: {len(pairs)} blind points   chi2 = {chi2:.2f}   "
      f"chi2/dof = {chi2/len(pairs):.2f}")
print(f"  worst point: {worst[1]} at {worst[0]:+.2f} sigma")

# ---------- 2. LOCALIZATION LAW: r(lambda, f1) everywhere per-config data exists ----------
print("\n" + "=" * 70)
print("2. LOCALIZATION LAW — r(lambda, f1), partial (plaquette removed)")
print("=" * 70)
zs = []
npos = ntot = 0
for f in sorted(glob.glob('results/results1[4-9]*.json')):
    R = json.load(open(f))
    for t, e in R['ensembles'].items():
        if 'cfg' not in e:
            continue
        lam = np.array(e['cfg']['lam']); f1 = np.array(e['cfg']['f1'])
        pl = np.array(e['cfg']['plaq']); n = e['n']
        r = partial(lam, f1, pl)
        s = z(r)*np.sqrt(n - 4)
        zs.append(s); ntot += 1; npos += (r > 0)
        print(f"  {f.split('/')[-1][:-5]}:{t:4} (b={e['beta']}, L={e['L']:2}, n={n:3})  "
              f"r = {r:+.3f}  ({s:+.1f} sigma)")
print(f"\n  sign-consistency: {npos}/{ntot} positive   "
      f"combined (Stouffer) = {sum(zs)/np.sqrt(len(zs)):+.1f} sigma")

# ---------- 3. THE LAMBDA LADDER: free + constrained RGZ fits ----------
print("\n" + "=" * 70)
print("3. LAMBDA_RGZ LADDER — 8-shell fits (free at G1; constrained elsewhere)")
print("=" * 70)


def shells_ph2(L):
    k = np.fft.fftfreq(L)*L
    kk = np.zeros((L,)*4); ph = np.zeros((L,)*4)
    for ax in range(4):
        s = [1]*4; s[ax] = L
        kk = kk + (k.reshape(s)**2); ph = ph + (2*np.sin(np.pi*k.reshape(s)/L))**2
    return np.array([ph[kk == s].mean() for s in range(1, 9)])


def load_pt(fname, tag):
    R = json.load(open(fname))
    e = R['ensembles'][tag]
    D = np.array(e['cfg']['D_shells'])
    return e, shells_ph2(e['L']), D.mean(0), D.std(0, ddof=1)/np.sqrt(e['n'])


def free_rgz(ph2, Dm, De):
    w = 1/De**2
    def c2(p):
        M1s, M2s, l4 = np.abs(p)
        sh = (ph2 + M1s)/(ph2**2 + M2s*ph2 + l4)
        Z = np.sum(w*Dm*sh)/np.sum(w*sh**2)
        return float(np.sum(w*(Dm - Z*sh)**2))
    best = None
    for x0 in ([0.6, 0.12, 0.04], [0.3, 0.3, 0.08], [1.0, 0.05, 0.02],
               [0.2, 0.1, 0.015], [0.4, 0.2, 0.1]):
        r = minimize(c2, x0, method='Nelder-Mead',
                     options={'xatol': 1e-7, 'fatol': 1e-9, 'maxiter': 6000})
        if best is None or r.fun < best.fun:
            best = r
    M1s, M2s, l4 = np.abs(best.x)
    return M1s, M2s, l4, best.fun


def constrained_lam(ph2, Dm, De, a, M1p, M2p, win=4.1):
    cut = ph2 <= win*a*a
    p, Dc, Ec = ph2[cut], Dm[cut], De[cut]
    w = 1/Ec**2
    M1l, M2l = M1p*a*a, M2p*a*a
    grid = np.linspace(0.0005, 0.4, 6000)
    cc = []
    for v in grid:
        sh = (p + M1l)/(p**2 + M2l*p + v)
        Z = np.sum(w*Dc*sh)/np.sum(w*sh**2)
        cc.append(float(np.sum(w*(Dc - Z*sh)**2)))
    cc = np.array(cc); i = int(cc.argmin()); ok = cc <= cc[i] + 1
    return (grid[i]**0.25/a, grid[ok].min()**0.25/a, grid[ok].max()**0.25/a,
            cc[i], int(cut.sum()))


# reference: G1 free fit
eG, phG, DmG, DeG = load_pt('results/results16.json', 'G1')
aG = ASQS[eG['beta']]
M1s, M2s, l4, c0 = free_rgz(phG, DmG, DeG)
M1p, M2p = M1s/aG**2, M2s/aG**2
print(f"  REFERENCE (free fit) G1 (b=2.4, L=20): chi2={c0:.2f}/4  "
      f"lambda = {l4**0.25/aG:.3f} sqs   M1p={M1p:.2f} M2p={M2p:.3f}")
print("  DISCLOSURE: all points below transfer M1,M2 from this single free fit;")
print("  the ladder is therefore consistency under a shared-parameter assumption,")
print("  not four independent measurements.")
for fname, tag in (('results/results19.json', 'J1'), ('results/results16.json', 'G2'),
                   ('results/results18.json', 'H1'), ('results/results18.json', 'H2')):
    e, ph2, Dm, De = load_pt(fname, tag)
    a = ASQS[e['beta']]
    lam, lo, hi, c2v, ns = constrained_lam(ph2, Dm, De, a, M1p, M2p)
    print(f"  {tag} (b={e['beta']}, L={e['L']}, a={a}): lambda = {lam:.3f} "
          f"[{lo:.3f}, {hi:.3f}]  (chi2={c2v:.2f}, {ns} shells)")
# free fit at J1 for the disclosure
eJ, phJ, DmJ, DeJ = load_pt('results/results19.json', 'J1')
M1j, M2j, l4j, c0j = free_rgz(phJ, DmJ, DeJ)
print(f"  free fit at J1 (disclosure): lambda = {l4j**0.25/ASQS[2.6]:.3f}, "
      f"chi2={c0j:.2f}/4, M1^2={M1j:.1f} (degenerate) — does NOT independently pin 1.0")

# ---------- 4. dead-conjecture count (this record) ----------
print("\n" + "=" * 70)
print("4. REFUTED IN THIS RECORD (runs 12-19)")
print("=" * 70)
for line in ("H-I: mass grows toward the horizon (F5, run 14)",
             "onset-at-death: mass appears when coherence dies (run 16/19, P17)",
             "energy-floor: horizon is the lowest-energy place (P11, run 15)",
             "E ~ 1/m per configuration (exploratory test, 2026-08-19)",
             "pure-GZ gamma + its a^2.4 'artifact law' (8-shell refit, run 16)",
             "IR-specificity of the horizon coupling (P2/F3, run 12)",
             "f3 and f3n non-spectral coherence gates (F2 run 11 / F4 run 12)"):
    print(f"  - {line}")
print("  (Part I's 13 refuted conjectures are separate; see companion record.)")

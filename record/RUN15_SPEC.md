# Run 15 — pre-registered specification: the INVERSE hypothesis (H-II)

*Registered 2026-08-19 08:34, before any run-15 data. Git-frozen locally at the
author's request; precedence timestamped by the Kaggle push. Fresh seeds; no
run-14 configuration is reused.*

## Provenance

Run 14 refuted H-I (mass grows toward the horizon): F5 fired — in the one
scaled ensemble (D1, dead vacuum) the FAR half carried the LARGER gamma
(0.499 [0.487,0.511] vs 0.451 [0.435,0.465]). The inverse was immediately proposed: 
MASS GROWS WITH DISTANCE FROM THE HORIZON — the Gribov interior is the
heavy region, the wall is light. Exploratory quartiles on the same data (labeled
exploratory, recorded in the analysis log) were NOT monotone (0.489, 0.412,
0.498, 0.505 near->far): the half-split signal may be real structure or n=16
noise. Run 15 decides on fresh data with doubled statistics.

## Design

Fresh-seed dead-vacuum ensembles (the only phase where a scale exists):
  E1 (2.2, L=10, n=128)  — the D1 replica, doubled statistics
  E2 (2.2, L=12, n=64)   — same coupling, bigger box (f ~ 0.38)
  E3 (2.3, L=14, n=48)   — different coupling, dead by grid (f ~ 0.38)
Per config: lambda_min, f1, D on shells k2=1..4, plaquette — per-config arrays
stored (run-14 instrumentation).

## Pre-registered predictions

P6 (primary, halves): in EVERY ensemble, the FAR half (lambda > median) fits a
    LARGER M2 scale than the NEAR half, with 1-sigma intervals separated.
P7 (secondary, structure): quartile gammas are monotone non-decreasing in lambda
    (Kendall tau > 0 across the 4 bins in at least 2 of 3 ensembles). If P6
    passes but P7 fails, the claim is "far-heavier" only, with unresolved fine
    structure — stated at that strength.
P8 (standing): blind r(k2=1) from the frozen activation curve for all three
    (dead phase: |r| < 0.25 expected).

## Falsifiers

F6: any ensemble with NEAR >= FAR (separated), or all three with overlapping
    intervals -> H-II dies; the run-14 half-split was noise; per-config
    lambda-dependence of the scale is closed in this programme (two named
    hypotheses, both refuted — no third variant will be registered).
F7: P6 passes at 2.2 but fails at 2.3 -> coupling-specific artifact, not a law.

## Verdict language (fixed now)

P6+P7 pass -> "within the dead-vacuum phase, the infrared scale grows with
distance from the Gribov horizon" — the per-configuration mass structure exists,
with the sign the run-14 data indicated.
F6 fires -> reported with full prominence; the mass remains a purely collective
(ensemble-level) phenomenon in this programme.

## ADDENDUM — registered 2026-08-19 08:48, kernel status RUNNING, no result seen

Two predictions from run-14 exploratory mining (labeled exploratory there),
frozen now so run 15 scores them blind:

P9  (the Z-law, from Hint 1): the per-config partial correlation of lambda_min
    with the BAND amplitude Z = mean over shells of D*phat^2 is stronger than
    with D(k2=1) alone. Run-14 values: ON plateau r_Z ~ -0.615, dead-phase
    r_Z = -0.24 (ambiguous: soft gate vs noise). Blind fork for run 15's three
    DEAD ensembles: (a) hard gate -> |r_Z| < 0.12; (b) soft gate -> r_Z in
    [-0.35, -0.15] in all three. Whichever branch the data takes, it decides
    the dead-phase Z-law; a mixture (some ensembles a, some b) refutes both.

P10 (wall-localizes, from Hint 2): r(lambda, f1) > 0 (partial, plaquette
    removed) in ALL THREE ensembles — configurations nearer the horizon carry a
    LESS coherent lowest mode. Run-14 saw +0.21..+0.50 in four of four
    ensembles. F10: any ensemble with r(lambda, f1) < 0 at >2 sigma refutes the
    per-config localization claim.

Also recorded: the trans-horizon count (configs with lambda < -10*tol) per
ensemble is reported as a new observable (no prediction; first measurement).

## ADDENDUM 2 — registered 2026-08-19 12:13, v2 kernel RUNNING, no v2 result seen
## (v1 crashed before writing output; its log contained no lambda-plaquette data)

P11 (energy-floor flip, from run-14 exploratory): in the DEAD phase the horizon
    is the energy floor — r(lambda, plaq) < 0 in all three ensembles (near-
    horizon configs have HIGHER plaquette = lower action). Run-14: ON ensembles
    showed the opposite sign (+0.16..+0.23), dead showed -0.16; the claim is the
    SIGN FLIP is a gate phenomenon. F11: any dead ensemble with r(lambda, plaq)
    > 0 at >2 sigma refutes the floor-migration picture.

# DIM-01 — pre-registered audit of the dimensional fixed point

*Registered 2026-08-24 18:11 UTC. The default flow had been run at this point
(6.121 -> 3.066 in 179 steps); the sweep in `dim_sweep.py` was launched before
this file was written and **its output had not been read** when the predictions
below were fixed. Nothing in the predictions is retro-fitted to the sweep. One
amendment was made after reading the first sweep, and it changed D7's verdict —
see Amendment 1, which is the first thing to read here.*

## Amendment 1 (2026-08-24 18:28 UTC) — the scoring statistic, changed after seeing data

The predictions below were registered against `d_inf`, an exponential-approach
fit. The first sweep exposed two defects in that fit, both in the analysis and
neither in the model:

1. a fixed 120-step tail returned *nothing* for any run shorter than 140 steps,
   which silently deleted three of the five floor-block runs — including both
   fast ones — from D7's own regression;
2. an unbounded tau grid fitted a large slow mode to a nearly flat tail and
   extrapolated below anything the run visited (`d_inf = 2.305` for a run that
   ended at 3.111).

The fit is repaired (adaptive tail, tau capped at the tail length) and demoted
to a cross-check. **The primary statistic is now `d_halt`: the measured
dimension at the step where the rate crossed the stopping tolerance, on runs
that actually converged (`rate_vanished`); runs that hit the step budget are
excluded and named.** All limits and all D7 branch readings are unchanged.

**This change was made after reading data, and it reversed D7.** On the broken
statistic D7's slope was 0.045 — which by the branch reading below would have
read as "the halt is NOT set by the floor", i.e. the emergence claim surviving.
On the repaired statistic it is 0.947, and the claim fails. Anyone auditing this
record should know that the one post-hoc methodology change here is the one that
made the headline result negative rather than positive, and can check it
directly: the broken regression ran on 3 of 5 points, the repaired one on the 4
that converged, and the raw halting values were the same in both sweeps
(2.103, 2.551, 3.043, 3.518 against floors 2.0, 2.5, 3.0, 3.5).

## What is being audited

A toy model in which dimension is not a parameter but a statistic: a cloud of
correlated events whose effective dimension is measured at every step by a
fixed-radius MLE, evolving under a rate built from two pieces of textbook
stability analysis —

- **where decay starts.** In d spatial dimensions the two-body effective
  potential is `-k/r^(d-2) + L^2/(2 mu r^2)`; circular orbits are stable iff
  d < 4, and the Coulomb problem falls to the centre for d >= 4. Structure
  above the threshold is erased. Entered as `d_collapse`.
- **where decay can stop.** The rate is multiplied by the propagating
  degree-of-freedom count of the interaction doing the erasing —
  `D(D-3)/2` in linearised GR, whose zero at D = 3 is a theorem rather than a
  dial. Entered as the `dof` floor, with a `linear` alternative `(d - d_grav)`
  whose zero can be moved on purpose.

Reported behaviour under test: a decay from ~6.1 that halts near 3, with the
instability rate switching itself off, and the claim that **nothing in the code
said "go to 3."**

## The claim as it can actually be tested

The strong reading ("3 is nowhere in the code") is false by inspection: the
`dof` floor vanishes at 3 by construction. The testable reading is

> the halting dimension is set by the floor and by nothing else — not by the
> soft parameter, not by the step size, not by the initial condition, not by
> the upper threshold — and the flow reaches it instead of stalling at the
> upper threshold or running through to zero.

If that holds, the model's content is a mechanism (a decay that self-arrests
at the zero of the graviton count) and the number 3 is inherited from
linearised GR, not produced by the dynamics. If it fails, the number is a
parameter artifact and there is no result.

## Pre-registered predictions

Fixed points were registered as `d_inf` from the exponential-approach fit over
the last 120 steps, on the reasoning that the rate vanishes linearly at the
floor so the run never literally arrives and a raw final value would quote the
step budget. Amendment 1 replaced this with `d_halt` on converged runs; read
`d_inf` below as `d_halt` throughout. The concern was real but the cure was
worse: `rate_vanished` fires at a rate 1600x below its starting value, which
pins the halt to within +0.054 of the floor — tighter than the fit's own scatter.

- **D1 (self-arrest).** In at least 4 of 5 seeds the default flow terminates on
  `rate_vanished`, not `max_steps`, and not by the estimator starving.
- **D2 (seed stability).** `d_inf` spread across 5 seeds < 0.15.
- **D3 (soft parameter).** `d_inf` spread across w in {0.25, 0.5, 0.75, 1.0}
  < 0.15. This is the model's one soft parameter; if the fixed point rides on
  it, the fixed point is not a result.
- **D4 (step size).** `d_inf` spread across eta0 in {0.3, 0.6, 1.2} < 0.15.
  A drift here means the halt is a discretisation artifact — the flow
  overshooting a floor it approaches in finite steps.
- **D5 (basin).** Starting from d0 in {4.6, 6.1, 8.3, 10.0} all give the same
  `d_inf` to within 0.20, and in particular a start at 8.3 or 10.0 does **not**
  stall at the upper threshold near 4. Registered separately: a start at 2.6,
  below the floor, does **not** rise — this model has no restoring term, so the
  fixed point is expected to be one-sided (an absorbing halt, not a two-sided
  attractor). Confirming that is a limitation, not a success.
- **D6 (upper threshold irrelevance).** With the `linear` floor held at
  d_grav = 3, `d_inf` spread across d_collapse in {3.5, 4.0, 4.5, 5.0} < 0.20.
  d_collapse should set how fast the decay runs, not where it ends.
- **D7 (the decisive control).** Regressing `d_inf` on d_grav over
  {1.5, 2.0, 2.5, 3.0, 3.5} with the `linear` floor gives a slope. Registered
  reading, fixed now:
  - slope > 0.8 -> **the halt is the input floor.** The number 3 is inherited
    from the graviton-count zero, and the emergent content of the model is the
    approach only. This is the outcome to be reported as the headline if it
    occurs, replacing the "chose three on its own" language wherever it
    appears.
  - slope < 0.3 -> the halt is set by something else in the dynamics and the
    emergence claim survives in a much stronger form than expected.
  - 0.3–0.8 -> mixed; the floor is one of several things setting the halt and
    no clean statement is made.

## Falsifiers

- **F-D1.** D3 or D4 failing (spread > 0.15) kills the fixed point outright:
  a halting value that rides on the soft parameter or the step size is a
  numerical artifact and the run is reported as such.
- **F-D2.** D1 failing in the default configuration — the flow running to zero,
  or stalling above 4 — means this reconstruction does not reproduce the
  reported behaviour, and nothing about the original model can be concluded
  from it either way.
- **F-D3.** If the calibration's max residual exceeds 0.10 over k = 1..8, the
  instrument is not linear over the range in use and no dimension quoted here
  is trustworthy to the precision claimed.

## Registered caveats

1. This is a reconstruction from a description, not the original code. It can
   only test whether the *mechanism as described* has the claimed property. If
   the original halts at 3 by some route other than a hand-placed zero, this
   audit does not touch it — but then D7 run on the original will show slope
   ~0, and that is the check to run there.
2. There is no time, no quantum mechanics, and no matter in this model. The
   stability analysis that supplies both thresholds is about atoms and orbits;
   nothing in the simulation contains either. The thresholds enter as numbers,
   with the physics standing outside the code as their justification.
3. The dimension measured is the dimension at one scale — the flow's own
   measurement ball. A state can read 3 there and read something else at every
   other scale, and the frozen state does exactly that. "Fixed point at 3"
   therefore means "fixed point of the dimension seen in this window", which is
   weaker than "a 3-dimensional space" and must not be quoted as the latter.
   The N = 30000 control (an exact d = 3 cloud, which reads ~3 at every window)
   is carried alongside every scale-resolved measurement so the difference
   stays visible.

---

## Scored 2026-08-24 18:30 UTC (appended after the sweep; the text above is unchanged)

D1 PASS (5/5 `rate_vanished`). D2 PASS (0.041). D3 PASS (0.038). D4 PASS
(0.027). D5 PASS (4 starts from 4.6 to 9.9 land in 3.035–3.111, none stalling
at the upper threshold); the below-floor start halts after one step, confirming
the halt is one-sided as registered. D6 PASS (0.049). F-D1, F-D2, F-D3 all
unfired.

**D7 FIRES in the first branch: slope 0.947, R² 0.9997.** Per the reading fixed
above, the halting dimension is the input floor, the "chose three on its own"
language is retired, and the model's content is the approach. Numbers and
discussion in `README.md`; raw output in `TALLIES.txt`.

## Grilled 2026-08-24 19:05 UTC — `grill.py`, adversarial review of this audit

Six objections run against DIM-01 itself; full numbers in `TALLIES.txt` section
4 and discussion in `README.md`.

- **G1 (attack failed).** The per-probe spread is real, not estimator noise:
  split-half reliability 0.77–0.96, and near the floor SD 0.519 = signal 0.510
  + noise 0.100.
- **G2/G3 (my own claim overturned).** D5's reading that "dispersion is what
  carries the flow past the upper threshold, and nothing puts it in by hand"
  is wrong as a causal claim. Deleting the spread changes nothing (6.121 ->
  3.075 in the same 179 steps); driving w to 0.030 changes nothing. Removing
  BOTH stalls the flow at 3.568. The pass-through is over-determined; the
  dispersion is sufficient but not necessary.
- **D7 and D2–D6 sharpened.** With Gamma = U x G and only G carrying a zero,
  the halt cannot depend on w, eta0, N or the seed, and must sit at G's zero.
  These passes are structural, readable off the rate without running anything.
  D7 is properly read as "nothing intervenes between the flow and the zero".
- **G4 (headline number corrected).** The halt is tolerance-limited at
  gamma_stop = 1e-3; it converges to 3.040 by 1e-4. The 3.066 quoted at
  registration is +0.026 high.
- **G5 (new, and stronger than D7).** The halt is a property of the measuring
  window, not of the geometry: windows R = 1.0, 1.4, 2.0 all halt at ~3.06 and
  all leave different frozen states (read at R = 1.0 they are 3.18, 4.17,
  4.75). "Froze at 3" is a statement about the instrument.
- **G6 (defect in this reconstruction).** The rate compares per-probe
  dimensions to d_collapse = 4 while every reported dimension is the pooled
  estimate; the Jensen gap between them is +0.194, so the upper threshold bites
  at pooled 3.81. Affects where the decay starts, not where it stops; no
  conclusion changes, but "passes through 4" should read "through 3.8".

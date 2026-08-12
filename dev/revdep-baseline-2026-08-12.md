# Downstream baseline before aligning `clarabel_control()` defaults

Recorded 2026-08-12, before any change to the defaults in `R/clarabel.R`.

The point of this file is to make a later claim falsifiable. If we ever flip the
chordal decomposition defaults (see "Deferred: chordal" below), the only honest
way to judge the result is to diff against a state that was written down *before*
the change, including the failures that were already there.

## Versions and platform

| item | value |
|:-----|:------|
| clarabel (this package) | 0.11.2 |
| Clarabel.rs (vendored) | 0.11.1, features `sdp-r` + defaults (no `faer-sparse`, no pardiso) |
| CVXR (installed) | 1.9.1.9017 |
| R | 4.6.0 (2026-04-24) |
| platform | aarch64-apple-darwin23, macOS Tahoe 26.5.1 |

## clarabel's own test suite

Full local suite, `tinytest::test_package("clarabel", at_home = TRUE)`:
**93 tests, 0 failures.**

| file | tests |
|:-----|------:|
| test_api_dimension.R | 9 |
| test_basic_eq_constrained.R | 4 |
| test_basic_expcone.R | 5 |
| test_basic_genpowcone.R | 2 |
| test_basic_lp.R | 13 |
| test_basic_powcone.R | 2 |
| test_basic_psd.R | 9 |
| test_basic_qp.R | 15 |
| test_basic_socp.R | 9 |
| test_basic_unconstrained.R | 3 |
| test_empty_conespec.R | 1 |
| test_mixed_conic.R | 2 |
| test_presolve.R | 7 |
| test_sdp_chordal.R | 12 |
| **Total** | **93** |

`test_basic_psd.R` and `test_sdp_chordal.R` (21 tests) are the ones exposed to a
chordal default change. The CRAN subset is only `test_basic_lp.R` +
`test_basic_qp.R` (28 tests); everything else is `at_home()`-gated, so a plain
`R CMD check` will *not* catch an SDP regression. Run with `at_home = TRUE`.

## Blast radius

A `clarabel_control()` default reaches CRAN through two hops.

**Hop 1 — direct reverse dependencies of clarabel (8):**
`aramappings`, `BsplineQuantReg`, `CVXR`, `diffcp`, `EmpiricalDynamics`,
`optweight`, `pensynth`, `ROI.plugin.clarabel`.

**Hop 2 — reverse dependencies of CVXR (36 checked by revdepcheck):** CVXR calls
bare `clarabel::clarabel_control()` and pins nothing, so whatever we ship becomes
CVXR's behavior on its next release.

Three packages (`aramappings`, `BsplineQuantReg`, `EmpiricalDynamics`) appear in
both hops. Net: **~40 distinct CRAN packages downstream.**

## revdepcheck baseline (CVXR)

From `~/research/cvxr/new_design/revdep_check_fresh.log`, run 2026-06-08,
`OK: 33, BROKEN: 3`, 36 packages, 12 min on 14 workers.

**Important caveat on attribution.** That run varied *CVXR* (1.8.2 → 1.9.1) with
clarabel held fixed. It is not a clarabel baseline. Its three "new problems" are
CVXR-attributable and must not be read as clarabel regressions. What it is good
for is the absolute per-package E/W/N counts below — the noise floor that any
future run has to be compared against.

| package | version | E | W | N | note |
|:--------|:--------|--:|--:|--:|:-----|
| aramappings | 0.2.0 | 1 (+1) | 0 | 0 | newly broken; cause is `parallel::makeCluster` port binding, environmental |
| ccar3 | 0.1.0 | 1 | 1 | 0 | pre-existing |
| cuadramelo | 1.0.0 | 2 | 1 | 1 | pre-existing |
| DebiasInfer | 0.2.1 | 0 | 0 | 0 | |
| DiSCos | 0.1.3 | 1 | 0 | 0 | failed to check — R CMD check timed out, in both |
| EmpiricalDynamics | 0.1.3 | 0 | 0 | 0 | also a direct clarabel revdep |
| EpiForsk | 0.2.2 | 1 | 0 | 0 | pre-existing |
| fairml | 0.9.1 | 0 | 0 | 0 | |
| fdaSP | 1.1.2 | 0 | 0 | 0 | |
| filling | 0.2.4 | 0 | 0 | 0 | |
| fungible | 2.4.7 | 0 | 0 | 0 | |
| graphicalExtremes | 0.3.5 | 1 | 0 | 0 | pre-existing |
| gtDesign | 0.1.0 | 0 | 0 | 0 | |
| hedgedrf | 1.0.1 | 1 | 1 | 0 | pre-existing |
| HonestDiD | 0.2.8 | 0 | 0 | 0 | |
| kantorovich | 3.2.0 | 3 | 1 | 1 | pre-existing |
| MaximinInfer | 2.1.0 | 0 | 0 | 0 | |
| mlr3fairness | 0.4.0 | 0 | 0 | 0 | |
| PlackettLuce | 0.4.5 | 0 | 0 | 0 | |
| PortfolioAnalytics | 2.1.2 | 0 | 0 | 0 | |
| portfolioBacktest | 0.4.2 | 0 | 0 | 0 | |
| pvEBayes | 0.2.2 | 0 | 0 | 0 | |
| rclsp | 2.0.0 | 0 (+1) | 0 | 0 | newly broken; package dependencies ERROR |
| Rdimtools | 1.1.4 | 0 | 0 | 0 | |
| RobustIV | 0.3.1 | 0 | 0 | 0 | |
| savvyGLM | 0.1.4 | 0 | 0 | 0 | |
| scpi | 4.0.0 | 0 | 0 (+1) | 0 | newly broken; installation WARNING |
| SIHR | 2.1.1 | 0 | 0 | 0 | |
| SLSEdesign | 0.1.0 | 0 | 0 | 0 | |
| spBPS | 2.0.1 | 0 | 0 | 0 | |
| spStack | 1.1.3 | 0 | 0 | 0 | |
| starnet | 1.0.2 | 1 | 0 | 0 | pre-existing |
| tramnet | 0.0.991 | 0 | 0 | 0 | |
| transreg | 1.0.7 | 1 | 0 | 0 | pre-existing |
| wdnet | 1.2.4 | 0 | 0 | 0 | |
| WRI | 0.2.3 | 0 | 0 | 0 | |

Artifacts live in `~/research/cvxr/new_design/CVXR/revdep/`
(`README.md`, `problems.md`, `failures.md`, `cran.md`, `data.sqlite`).

## Changes being made now

These are believed inert or strictly corrective. Recorded here so that if
something downstream *does* move, we know the candidate causes.

1. `direct_solve_method`: `"qdldl"` → `"auto"`, options narrowed to
   `c("auto", "qdldl")`. Provably a no-op in our build — without `faer-sparse`,
   `AutoDirectLDLSolver::new` reduces to `QDLDLDirectLDLSolver::new(KKT, Dsigns,
   settings, None)`, the identical call the `"qdldl"` branch makes. Also removes
   `"cholmod"` (never existed in Clarabel.rs; it is a Clarabel.jl option) and
   `"mkl"` (needs the `pardiso-mkl` feature we do not build; now hard-errors
   under v0.11's construction-time `validate()`).
2. `reduced_tol_infeas_abs`: `5e-5` → `5e-12`. Stale value. Upstream changed this
   in `5e94148` (2025-01-07, first released in Clarabel.rs v0.10.0); we were
   never resynced. Per upstream's own note, smaller means *less* margin required,
   so we have been stricter than upstream. Aligning makes
   `AlmostPrimalInfeasible` / `AlmostDualInfeasible` marginally easier to reach,
   which CVXR maps to `INFEASIBLE_INACCURATE` / `UNBOUNDED_INACCURATE`. Affects
   borderline problems only.
3. Documentation typos, no code change: `tol_ktratio` (doc said `1e-7`, code and
   upstream are `1e-6`), `iterative_refinement_reltol` (doc said `1e-12`, code
   and upstream are `1e-13`).
4. A drift test asserting `clarabel_control()` against
   `DefaultSettings::<f64>::default()` field by field, with the chordal
   deviations allowlisted and commented.

## Chordal decomposition: measured, then aligned with upstream

The first draft of this document argued for keeping the four
`chordal_decomposition_*` settings off, against upstream, on the grounds that
enabling them changes reported PSD duals for ~40 downstream packages while
buying the warm-start work nothing. That argument was made without measuring the
speedup. The measurement killed it.

Benchmark: linear matrix inequality `F0 + sum_i x_i Fi >= 0`, all `Fi` banded
with half-bandwidth 2 (so the aggregate pattern is banded, hence chordal), 5
variables. Script in this repo's scratchpad; regenerate with `bench_chordal.R`.

Re-run with the dense construction the vignette uses (`dev/bench_chordal.R`
was the sparse first pass; `sweep.R` output below is what Section 3.2 reports):

| n | PSD cone length | chordal off | chordal on | speedup | status off / on | obj abs diff | max abs dz |
|---:|---:|---:|---:|---:|:---|---:|---:|
| 20 | 210 | 0.013s | 0.001s | 13x | Solved / Solved | 1.8e-09 | 8.9e-06 |
| 40 | 820 | 0.420s | 0.002s | 210x | AlmostSolved / Solved | 1.4e-07 | 6.7e-01 |
| 60 | 1830 | 4.393s | 0.004s | 1098x | Solved / Solved | 7.2e-08 | 6.3e-04 |
| 80 | 3240 | 29.855s | 0.004s | 7464x | Solved / Solved | 4.0e-08 | 2.3e-05 |
| 120 | 7260 | 424.978s | 0.008s | 53122x | Solved / AlmostSolved | 1.2e-07 | 8.3e-04 |

`length(z)` is identical either way (1830 at n=60), confirming that
`decomp_reverse` reassembles at the original dimensions.

Note the reduced-accuracy terminations, which fall on the undecomposed side
at n=40 and on the decomposed side at n=120. Decomposition produces a
genuinely different numerical problem, so neither path dominates the other on
accuracy; a seed scan at n=40 showed the split going both ways across
instances. This is reported as-is in the vignette rather than tuned away by
choosing a well-behaved seed — it is the reason to check `status` instead of
assuming `Solved`.

The "off" column grows roughly like n^6; the "on" column is close to linear.
This is not a tuning knob, it is a change of asymptotic complexity. Defaulting it
off is not the conservative choice — it is a trap for anyone solving a sparse
semidefinite program who does not know the flag exists. All four settings are
now aligned with upstream.

What the earlier argument got right, and what it got wrong:

- Right: PSD **dual values** genuinely move (up to 0.67 at n=40). That is the
  real risk and it remains the thing to watch downstream.
- Wrong: the claim that the flip buys warm starts nothing, and therefore costs
  nothing to refuse. Chordal decomposition only engages when PSD cones are
  present, and updates are blocked by `is_chordal_decomposed()` — whether a
  decomposition *actually happened*, not whether it was enabled. Linear,
  quadratic and second-order cone problems, which are the bulk of CVXR's
  warm-start traffic, are unaffected and stay updatable. The two features only
  collide on SDPs, where the measurement says decomposition is worth orders of
  magnitude and warm starting a small constant.

`clarabel_control()` now matches `DefaultSettings::<f64>::default()` on all 43
fields, asserted field by field in `inst/tinytest/test_control_defaults.R`. The
only remaining departure is the *set of accepted values* for
`direct_solve_method`, which is constrained by the Cargo features this package
builds with, not a choice.

## BLAS: reference vs vecLib

Semidefinite support is built with the `sdp-r` feature, so `clarabel.so` links
`libRblas`/`libRlapack` and the dense linear algebra comes from R. Ran the
Section 3.2 sweep under both macOS BLAS libraries by swapping the
`libRblas.dylib` symlink (`dev/blas_compare.R`, driven by a script that
restores the reference BLAS on any exit path).

| n | reference off | vecLib off | reference on | vecLib on |
|---:|---:|---:|---:|---:|
| 20 | 0.015s | 0.016s | 0.001s | 0.002s |
| 40 | 0.429s | 0.444s | 0.003s | 0.003s |
| 60 | 4.475s | 4.576s | 0.003s | 0.004s |
| 80 | 30.078s | 30.535s | 0.005s | 0.005s |
| 120 | 426.752s | 429.727s | 0.008s | 0.009s |

**No speed difference** — every ratio is 0.94x to 1.00x, i.e. vecLib is
marginally slower, within run-to-run noise. This contradicts what the vignette
originally asserted, which was that an optimized BLAS would shrink the chordal
speedup ratios. That was written as a plausible inference without measuring it.

The swap was verified to take effect: `extSoftVersion()` reported
Accelerate's `libBLAS.dylib` during the second run, `libRlapack.dylib` links
`libRblas.dylib` dynamically, and `clarabel.so` has undefined symbols for
`dgemm`, `dsyrk`, `dpotrf` and `dsyevr`. The dense work is on the n x n cone;
the cost is the KKT factorization, which scales with the vectorized length
n(n+1)/2 (7260 rows at n=120).

**Results do differ.** Objectives agree to ~1e-8, but duals move by up to
1.27e-01 (n=40, decomposed), and the n=40 decomposed solve changed status from
`Solved` under the reference BLAS to `AlmostSolved` under vecLib. Iteration
counts were identical except n=60 decomposed (16 vs 15).

Relevant to the CRAN thread-count question below: Accelerate threads via Grand
Central Dispatch, so a user on vecLib could exceed two cores where the
reference BLAS would not.

## Still to verify

1. CVXR's `test-clarabel-warm-start.R` and SDP-touching tests against this build.
2. revdepcheck for CVXR 1.9.2 (in progress separately) — compare against the
   per-package table above, remembering that its three "new problems" predate
   this work.
3. Of the ~40 downstream packages, which actually construct PSD constraints.
   That is the true exposed set for the dual-value change, and it is much
   smaller than 40.
4. `max_threads = 0` lets the solver choose its own thread count. Worth checking
   against CRAN's two-core policy, since this build reaches BLAS/LAPACK through
   R and QDLDL is single-threaded — likely inert, but unverified.

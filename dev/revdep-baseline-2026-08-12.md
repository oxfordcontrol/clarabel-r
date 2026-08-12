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

## Deferred: chordal decomposition

We are keeping `chordal_decomposition_enable = FALSE`,
`chordal_decomposition_merge_method = "none"`, `chordal_decomposition_compact =
FALSE`, `chordal_decomposition_complete_dual = FALSE`, against upstream, which
has defaulted all four to on / `"clique_graph"` since v0.9.0.

Reasons, so a future reader does not have to reconstruct them:

- The change is not observable as an error. `decomp_reverse` reassembles at
  `self.init_dims`, so `z` and `s` come back at original dimension and ordering
  and CVXR's offset arithmetic in `clarabel_extract_dual_value` stays valid. What
  changes is the *values*: `psd_completion()` runs, and a PSD completion is not
  unique. Downstream sees different numbers in `dual_value` with no warning.
- ~40 CRAN packages are downstream, several with stored test expectations.
- CVXR's warm-start path already forces `chordal_decomposition_enable <- FALSE`
  itself, so the flip would buy the warm-start work nothing.
- `FALSE` is the conservative setting: it never silently changes anyone's duals
  and never blocks data updates.

This is a pinned decision, not drift. The drift test will fail on every future
Clarabel bump and put the question back in front of us.

**To revisit it, measure rather than argue.** The experiment:

1. Build clarabel with the four chordal defaults flipped.
2. Run this package's suite with `at_home = TRUE`; watch `test_basic_psd.R` and
   `test_sdp_chordal.R` (21 tests).
3. Solve a chordally sparse SDP through `clarabel()` both ways and diff `x`,
   `z`, `s`, `obj_val`, and iteration count. Expect `x` and `obj_val` to agree to
   tolerance and the PSD block of `z` to differ.
4. Run CVXR's `test-clarabel-warm-start.R` and its SDP-touching tests against
   both builds.
5. Of the ~40 downstream packages, identify which actually construct PSD
   constraints — that is the true exposed set, and it is much smaller than 40.

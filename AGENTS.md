# AGENTS.md — investr

@../CLAUDE.md

R package for **inverse estimation** (the calibration problem): given an
observed value of the response, estimate the corresponding unknown predictor
value. Exports: `invest()` (numerical/bootstrap-based, supports `lm`, `glm`,
`nls`, `lme`), `calibrate()` (closed-form solution for simple linear
regression only), `predFit()` (predictions/standard errors/bands, used
internally by `plotFit()`), and `plotFit()` (base-graphics plotting of a
fitted model with confidence/prediction bands).

investr is a long-lived CRAN package (on CRAN continuously since 2014) with
a companion R Journal article (Greenwell & Schubert Kabban, 2014). **Keeping
CRAN availability is the top priority** — unlike this author's other
packages (pdp, fastshap, vip, sure), investr is not migrating off CRAN.
r-universe is a secondary distribution channel, not a replacement.

## Branches & releases

- **`devel`** (default): all development and PRs. NEWS.md starts with
  `# investr (development version)`.
- **`main`**: stable, CRAN-released version only, tagged `vX.Y.Z`.
  r-universe (pinned to main in bgreenwell/bgreenwell.r-universe.dev) and the
  pkgdown site both build from main.
- Release: merge devel → main (`--no-ff`), tag, submit to CRAN, `gh release
  create`; then merge main back into devel and bump the devel NEWS heading.
- Shared fixes main needs immediately: commit to **main first, then merge
  main → devel**. Never cherry-pick devel → main.

## Dependency philosophy

Imports are minimal: `graphics, stats, utils`. `nlme` is a Suggests-only
dependency used solely by the `"lme"`-method code paths
(`invest.lme()`, `predFit.lme()`, and the internal `computeWaldInterval.lme`
/ `varY` / `makeZ` helpers in `R/utils.R`); both public entry points guard
with `requireNamespace("nlme", quietly = TRUE)` before dispatching, since a
caller can only have an `"lme"` object if nlme was already loaded to fit it.

## Commands

```bash
Rscript -e 'devtools::document()'                                            # after roxygen edits
Rscript -e 'pkgload::load_all("."); tinytest::run_test_dir("inst/tinytest")'  # full test suite
Rscript -e 'devtools::check(args = c("--as-cran"))'                          # before pushing
```

Tests use **tinytest** (not testthat) in `inst/tinytest/`. Internal
(non-exported) functions must be called with `investr:::` — unlike
testthat's `test_check()`, tinytest test scripts don't get automatic
insider access to the package namespace. Always add a NEWS.md entry; never
edit `man/` by hand.

## Architecture (R/)

- `invest.R` — `invest()` generic and its `lm`/`glm`/`nls`/`lme` methods:
  extract data, compute a point estimate by inverting the fitted model
  (`computeInverseEstimate.R`), then dispatch to the requested interval type.
- `calibrate.R` — `calibrate()`: closed-form calibration/inversion interval
  for the simple linear regression model only (no numerical search). Reads
  x/y directly from `stats::model.frame()`, not `model.matrix()` — the
  latter expands terms like `poly(x, degree)` into multiple orthogonal
  columns, which silently corrupted results before this was fixed (#48/#49).
- `computeInverseEstimate.R` — root-finding (`stats::uniroot()`) for the
  point estimate, per model class.
- `computeInversionInterval.R` — inversion (Fieller-type) confidence
  intervals, per model class.
- `computeWaldInterval.R` — Wald-type intervals via the delta method
  (`stats::numericDeriv()`), per model class.
- `computeBootReps.R` — parametric/nonparametric bootstrap for
  `interval = "percentile"`; calls `stats::simulate()`, which dispatches to
  the internal `simulate.nls()` (`R/simulate.nls.R`) when called from within
  the package's own namespace (it's intentionally unexported — see Gotchas).
- `predFit.R` — `predFit()` generic and its `lm`/`glm`/`nls`/`lme` methods:
  fitted values, standard errors, and confidence/prediction intervals
  (including Bonferroni/Scheffe adjustment for simultaneous intervals).
- `plotFit.R` — base-graphics plotting built on top of `predFit()`.
- `utils.R` — shared helpers (`makeData`, `makeX`, `makeZ`, `varY`,
  `sigma.lme`), mostly supporting the `"lme"` methods.

## Gotchas

- **`simulate.nls()` is intentionally unexported** (`@keywords internal`,
  no `@export`, no S3 registration). It's called internally via
  `stats::simulate()` from `computeBootReps.R`, which works because S3
  dispatch for unexported methods succeeds when the call originates from
  code *inside* the defining package's namespace — but fails for any
  external caller (including test scripts). tinytest tests must call
  `investr:::simulate.nls()` directly rather than going through the generic.
- **`calibrate.lm()` and `poly()`**: see the `calibrate.R` note above. Any
  future change to how `x`/`y` are extracted must keep
  `inst/tinytest/test_calibrate.R`'s `mlr.poly` case green.
- Rd `@seealso`/cross-reference tags must only point at *declared*
  (Imports/Suggests) packages — CRAN flags undeclared Rd xrefs as a NOTE.
  `nlstools::plotfit` and `lme4::bootMer` are referenced in prose (not as
  `[pkg::fun()]` links) for exactly this reason.

## CI / site

GitHub Actions (r-lib/actions v2): R-CMD-check + test-coverage on pushes/PRs
to main and devel; pkgdown builds from **main only** and deploys to
`gh-pages` (served at https://bgreenwell.github.io/investr/). The vignette
(`vignettes/introduction.Rnw`) is a legacy Sweave/LaTeX document tied to the
2014 R Journal article — it is intentionally not being converted to Rmd.

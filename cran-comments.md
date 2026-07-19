## Summary

This release is primarily bug fixes (several correctness issues in `invest()`,
`calibrate()`, and `predFit()` reported on GitHub, some dating back several
years), plus two small new features (`invest()` methods for Kaplan-Meier
survival curves, and bootstrap results that are now also valid `"boot"`
objects). See NEWS.md for the full list.

One fix is a behavior change worth flagging explicitly: `calibrate.lm()`
previously read the predictor from `model.matrix()`, which silently expands
terms like `poly(x, degree)` into orthogonal-polynomial columns and produced
incorrect estimates for such models without any warning. It now reads from
the model frame directly and errors informatively for predictor
transformations it doesn't support, rather than silently returning wrong
numbers. A reverse dependency check did not turn up any packages relying on
the old (incorrect) behavior.

## Test environments

* local macOS install, R 4.5.2
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release,
  oldrel-1)
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 2 notes

* "checking for future file timestamps ... NOTE / unable to verify current
  time" -- environmental, unrelated to the package.
* "checking S3 generic/method consistency ... NOTE" -- internal helper
  functions (never exported, not part of the public API) follow a
  `generic.class` naming convention for dispatch but aren't registered via
  NAMESPACE, since they're never called from outside the package. This NOTE
  has been present on prior CRAN releases of investr.

## Reverse dependencies

We checked 5 reverse dependencies (chemCal, envalysis, enveomics.R,
ggtrendline, GrowthCurveME), comparing R CMD check results across CRAN and
this release.

* We saw 0 new problems
* We failed to check 0 packages

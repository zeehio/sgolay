# sgolay 1.0.4.9000

# sgolay 1.0.4 (2026-07-09)

- Fix a crash on R >= 4.3 in `sgolayfilt()` when a plain (unclassed) filter
  coefficient matrix is passed as `p` (`dim(p) > 1` errored inside `&&` for
  any matrix; replaced with `is.matrix(p)`).
- `sgolayfilt()` now validates that a custom filter coefficient matrix `p`
  has an odd number of rows, instead of silently producing incorrect
  boundary values.
- `sgolayfilt()` now errors with a clear message when `x` is shorter than
  the filter length, instead of an obscure "subscript out of bounds" error.
- Added regression tests for the "auto" engine's filter-length threshold and
  the fft-to-filter NA fallback.

# sgolay 1.0.3 (2023-03-30)

- Fix warning on bitwise comparison in my_isok() macro
- Fix lossy conversion on indices from long to int

# sgolay 1.0.2 (2022-09-29)

Fix BugReports URL (http->https)

# sgolay 1.0.1

Updated DESCRIPTION file to address CRAN reviewers' comments.

# sgolay 1.0.0

Initial release

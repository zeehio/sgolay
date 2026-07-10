# New release 1.0.4

This release fixes a crash in `sgolayfilt()` on R >= 4.3 when a plain
(unclassed) filter coefficient matrix is passed as `p`: `dim(p) > 1`
evaluated to a vector of length 2 for any matrix, which errors inside `&&`
on R >= 4.3 ("'length = 2' in coercion to 'logical(1)'"). It also adds
input validation (odd-length filter matrix, `x` at least as long as the
filter) with clear error messages instead of obscure indexing failures,
and adds regression tests for previously untested code paths.

## Test environments

* local Ubuntu 24.04 install, R 4.3.3 (`R CMD check --as-cran`)
* GitHub Actions R-CMD-check, all passing:
  * macOS-latest, R release
  * windows-latest, R release
  * ubuntu-latest, R devel
  * ubuntu-latest, R release
  * ubuntu-latest, R oldrel-1

## R CMD check results

0 errors, 0 warnings, 0 notes on all of the above.

## Downstream dependencies

There are no downstream dependencies.

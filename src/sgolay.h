#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP cconvolve_circular(SEXP sx, SEXP sy, SEXP sconj);
SEXP cfilter(SEXP sx, SEXP sfilter);

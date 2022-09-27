#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP pad_and_convolve(SEXP sx, SEXP sy, SEXP sconj);
SEXP filter(SEXP sx, SEXP sfilter);

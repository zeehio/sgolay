#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP pad_and_convolve(SEXP sx, SEXP sy, SEXP sconj);
SEXP filter(SEXP sx, SEXP sfilter);

SEXP c_hclust(SEXP sn, SEXP sd, SEXP smethod, SEXP smembers);
SEXP c_hcass(SEXP sn, SEXP ia, SEXP ib);

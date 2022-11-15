#include "sgolay.h"

#include "hclust.h"

#define R_NO_REMAP
#include <Rconfig.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

SEXP c_hclust(SEXP sn, SEXP sd, SEXP smethod, SEXP smembers) {
  R_xlen_t n = INTEGER(sn)[0];
  SEXP sia = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP sib = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP scrit = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP sdisscpy = PROTECT(Rf_duplicate(sd));
  SEXP sout = PROTECT(Rf_allocVector(VECSXP, 3));

  hclust(
    n,
    INTEGER(smethod)[0],
    INTEGER(sia),
    INTEGER(sib),
    REAL(scrit),
    REAL(smembers),
    (idx_t*) R_alloc(n, sizeof(idx_t)),
    (double*) R_alloc(n, sizeof(double)),
    REAL(sdisscpy),
    (int *) R_alloc(n, sizeof(int))
  );
  SET_VECTOR_ELT(sout, 0, sia);
  SET_VECTOR_ELT(sout, 1, sib);
  SET_VECTOR_ELT(sout, 2, scrit);
  UNPROTECT(5);
  return sout;
}

SEXP c_hcass(SEXP sn, SEXP ia, SEXP ib) {
  R_xlen_t n = INTEGER(sn)[0];
  SEXP siia = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP siib = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP sorder = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP sout = PROTECT(Rf_allocVector(VECSXP, 3));
  hcass2(
    n,
    INTEGER(ia),
    INTEGER(ib),
    INTEGER(sorder),
    INTEGER(siia),
    INTEGER(siib)
  );

  SET_VECTOR_ELT(sout, 0, sorder);
  SET_VECTOR_ELT(sout, 1, siia);
  SET_VECTOR_ELT(sout, 2, siib);
  UNPROTECT(4);
  return sout;

}


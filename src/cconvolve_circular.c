#include <string.h>
#include <stddef.h>

#define R_NO_REMAP
#include <Rconfig.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

#include "sgolay.h"

void fft_factor(int n, int *pmaxf, int *pmaxp);
Rboolean fft_work(double *a, double *b, int nseg, int n, int nspn, int isn,
                  double *work, int *iwork);

/* given sx as a real matrix and sy as a real vector
 * this function returns a real matrix `ans` of `nrow(sx) + length(sy) -1` rows and
 * ncol(sx) columns where:
 * ans[,j] = stats::convolve(sx[,j], sy, conj = sconj, type = "circular").
 *
 * Long dimensions are not supported
 */
SEXP cconvolve_circular(SEXP sx, SEXP sy, SEXP sconj)
{
  int conj, x_rows, x_cols, y_len, ans_rows, ans_cols;
  int i, j;
  double *xpad_r, *xpad_i;
  double *ypad_r, *ypad_i;

  int maxf, maxp;
  double *work;
  int *iwork;
  size_t smaxf;
  const size_t maxsize = ((size_t) -1) / 4;
  double *x, *y, *ans, *xvec;
  SEXP d;
  int num_elem;
  double scale_factor, dtmp;
  int increment;

  conj = Rf_asLogical(sconj);
  if (conj == NA_INTEGER || conj == 0) {
    conj = 0;
  } else {
    conj = 1;
  }


  switch (TYPEOF(sx)) {
  case INTSXP:
  case LGLSXP:
  case REALSXP:
    sx = Rf_coerceVector(sx, REALSXP);
    break;
  default:
    Rf_error("x is not real. It should be a real matrix");
  }
  PROTECT(sx);

  switch (TYPEOF(sy)) {
  case INTSXP:
  case LGLSXP:
  case REALSXP:
    sy = Rf_coerceVector(sy, REALSXP);
    break;
  default:
    Rf_error("y is not real. It should be a real vector");
  }
  PROTECT(sy);


  d = Rf_getAttrib(sx, R_DimSymbol);
  if (d == R_NilValue || Rf_length(d) > 2) {
    Rf_error("x should be a matrix");
  }

  x_rows = Rf_nrows(sx);
  x_cols = Rf_ncols(sx);
  x = REAL(sx);

  y_len = Rf_length(sy);
  y = REAL(sy);

  ans_rows = x_rows + y_len - 1;
  ans_cols = x_cols;

  /* prep work for fft */
  fft_factor(ans_rows, &maxf, &maxp);
  if (maxf == 0)
    Rf_error("fft factorization error");
  smaxf = maxf;
  if (smaxf > maxsize)
    Rf_error("fft too large");
  work = (double*)R_alloc(4 * smaxf, sizeof(double));
  iwork = (int*)R_alloc(maxp, sizeof(int));

  SEXP sans = Rf_allocMatrix(REALSXP, ans_rows, ans_cols);
  PROTECT(sans);
  ans = REAL(sans);

  /* prepare y for convolution: */
  ypad_r = (double*) R_alloc(ans_rows, sizeof(double));
  ypad_i = (double*) R_alloc(ans_rows, sizeof(double));

  /* set ypad to zero on the imag part */
  memset(ypad_i, 0, ans_rows*sizeof(double));
  /* set ypad to sy at the beginning of the real part */
  memcpy(ypad_r, y, y_len*sizeof(double));
  /* set ypad to zero on rest of the real part */
  memset(ypad_r + y_len, 0, (ans_rows-y_len)*sizeof(double));

  /* fft of ypad in place */
  fft_factor(ans_rows, &maxf, &maxp);
  fft_work(ypad_r, ypad_i, 1, ans_rows, 1, 1, work, iwork);

  /* conjugate ypad */
  if (conj) {
    num_elem = ans_rows;
    scale_factor = -1;
    increment = 1;
    F77_CALL(dscal)(&num_elem, &scale_factor, ypad_i, &increment);
  }

  xpad_i = (double*)R_alloc(ans_rows, sizeof(double));
  for (j=0;j<ans_cols;j++) {
    xvec = x + j*x_rows;
    xpad_r = ans + j*ans_rows;
    /* set sxpad_imag to zero */
    memset(xpad_i, 0, ans_rows*sizeof(double));
    /* set sans[first_elems,j] to zero */
    memset(xpad_r, 0, (y_len - 1)*sizeof(double));
    /* set sans[last_elems,j] to x[,j] */
    memcpy(xpad_r + y_len - 1, xvec, x_rows*sizeof(double));
    /* fft of sans */
    fft_factor(ans_rows, &maxf, &maxp);
    fft_work(xpad_r, xpad_i, 1, ans_rows, 1, 1, work, iwork);
    /* product of fft of sans and sypad, saving to fft of sans */

    /* There's no blas for complex-complex element wise product */
    /* xpad*ypad: (a + ib) (c + id) = (ac - bd) + i(ad + bc). */
    for (i = 0; i< ans_rows; i++) {
      dtmp = xpad_r[i]*ypad_r[i] - xpad_i[i]*ypad_i[i];
      xpad_i[i] = xpad_r[i]*ypad_i[i] + xpad_i[i]*ypad_r[i];
      xpad_r[i] = dtmp;
    }

    /* inverse of fft of sans */
    fft_factor(ans_rows, &maxf, &maxp);
    fft_work(xpad_r, xpad_i, 1, ans_rows, 1, -1, work, iwork);
  }
  /* Divide by the length */
  num_elem = ans_rows*ans_cols;
  scale_factor = 1./ans_rows;
  increment = 1;
  F77_CALL(dscal)(&num_elem, &scale_factor, ans, &increment);
  UNPROTECT(3);
  return(sans);
  }



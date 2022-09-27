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

/*
 * Long dimensions are not supported
 */
SEXP pad_and_convolve(SEXP sx, SEXP sy, SEXP sconj)
{
  int conj;
  int x_is_matrix, x_num_signals, x_length, x_numdim;
  int y_is_matrix, y_num_signals, y_length, y_numdim;
  int ans_num_signals, ans_length;
  int i, j;
  double *xpad_r, *xpad_i;
  double *ypad_r, *ypad_i;
  double *anspad_r, *anspad_i;

  int maxf, maxp;
  double *work;
  int *iwork;
  size_t smaxf;
  const size_t maxsize = ((size_t) -1) / 4;
  double *x, *y, *ans;
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
    Rf_error("x is not real. It should be a real matrix or vector");
  }
  PROTECT(sx);

  switch (TYPEOF(sy)) {
  case INTSXP:
  case LGLSXP:
  case REALSXP:
    sy = Rf_coerceVector(sy, REALSXP);
    break;
  default:
    Rf_error("y is not real. It should be a real vector or matrix");
  }
  PROTECT(sy);

  /* options:
   * (a) x is a vector, y is a vector
   * (b) x is a matrix, y is a vector
   * (c) x is a vector, y is a matrix
   * (d) x is a matrix, y is a matrix: Not supported
   */

  d = Rf_getAttrib(sx, R_DimSymbol);
  if (d == R_NilValue) {
    x_is_matrix = 0;
    x_num_signals = 1;
    x_length = Rf_length(sx);
  } else {
    x_numdim = Rf_length(d);
    if (x_numdim == 1) {
      x_is_matrix = 0;
      x_num_signals = 1;
      x_length = Rf_length(sx);
    } else if (x_numdim == 2) {
      x_is_matrix = 1;
      x_num_signals = Rf_ncols(sx);
      x_length = Rf_nrows(sx);
    } else if (x_numdim > 2) {
      Rf_error("x should be a vector or matrix");
    }
  }
  x = REAL(sx);

  d = Rf_getAttrib(sy, R_DimSymbol);
  if (d == R_NilValue) {
    y_is_matrix = 0;
    y_num_signals = 1;
    y_length = Rf_length(sy);
  } else {
    y_numdim = Rf_length(d);
    if (y_numdim == 1) {
      y_is_matrix = 0;
      y_num_signals = 1;
      y_length = Rf_length(sy);
    } else if (y_numdim == 2) {
      y_is_matrix = 1;
      y_num_signals = Rf_ncols(sy);
      y_length = Rf_nrows(sy);
    } else if (y_numdim > 2) {
      Rf_error("y should be a vector or matrix");
    }
  }
  y = REAL(sy);

  if (x_num_signals > 1 && y_num_signals > 1) {
    Rf_error("convolution not vectorized on both x and y, one of those should be a vector");
  }

  ans_length = x_length + y_length - 1;
  ans_num_signals = x_num_signals == 1 ? y_num_signals : x_num_signals;

  /* prep work for fft */
  fft_factor(ans_length, &maxf, &maxp);
  if (maxf == 0)
    Rf_error("fft factorization error");
  smaxf = maxf;
  if (smaxf > maxsize)
    Rf_error("fft too large");
  work = (double*)R_alloc(4 * smaxf, sizeof(double));
  iwork = (int*)R_alloc(maxp, sizeof(int));

  SEXP sans = Rf_allocMatrix(REALSXP, ans_length, ans_num_signals);
  PROTECT(sans);
  ans = REAL(sans);

  /* prepare the vector signal for convolution */
  if (y_num_signals == 1) {
    /* prepare y for convolution: */
    ypad_r = (double*) R_alloc(ans_length, sizeof(double));
    ypad_i = (double*) R_alloc(ans_length, sizeof(double));

    /* set ypad to zero on the imag part */
    memset(ypad_i, 0, ans_length*sizeof(double));
    /* set ypad to sy at the beginning of the real part */
    memcpy(ypad_r, y, y_length*sizeof(double));
    /* set ypad to zero on rest of the real part */
    memset(ypad_r + y_length, 0, (ans_length-y_length)*sizeof(double));
    /* fft of ypad in place */
    fft_factor(ans_length, &maxf, &maxp);
    fft_work(ypad_r, ypad_i, 1, ans_length, 1, 1, work, iwork);
    /* conjugate ypad */
    if (conj) {
      num_elem = ans_length;
      scale_factor = -1;
      increment = 1;
      F77_CALL(dscal)(&num_elem, &scale_factor, ypad_i, &increment);
    }
  } else {
    /* prepare x for convolution */
    xpad_r = (double*) R_alloc(ans_length, sizeof(double));
    xpad_i = (double*) R_alloc(ans_length, sizeof(double));
    memset(xpad_i, 0, ans_length*sizeof(double));
    /* set sans[first_elems,j] to zero */
    memset(xpad_r, 0, (y_length - 1)*sizeof(double));
    /* set sans[last_elems,j] to x[,j] */
    memcpy(xpad_r + y_length - 1, x, x_length*sizeof(double));
    /* fft of xpad in place */
    fft_factor(ans_length, &maxf, &maxp);
    fft_work(xpad_r, xpad_i, 1, ans_length, 1, 1, work, iwork);
  }
  double *helper_i = (double*)R_alloc(ans_length, sizeof(double));
  for (j=0;j<ans_num_signals;j++) {
    if (y_num_signals == 1) {
      xpad_i = helper_i;
      double *xvec = x + j*x_length;
      xpad_r = ans + j*ans_length;
      /* set sxpad_imag to zero */
      memset(xpad_i, 0, ans_length*sizeof(double));
      /* set sans[first_elems,j] to zero */
      memset(xpad_r, 0, (y_length - 1)*sizeof(double));
      /* set sans[last_elems,j] to x[,j] */
      memcpy(xpad_r + y_length - 1, xvec, x_length*sizeof(double));
      /* fft of sans */
      fft_factor(ans_length, &maxf, &maxp);
      fft_work(xpad_r, xpad_i, 1, ans_length, 1, 1, work, iwork);
      /* product of fft of sans and sypad, saving to fft of sans */

      /* There's no blas for complex-complex element wise product */
      /* xpad*ypad: (a + ib) (c + id) = (ac - bd) + i(ad + bc). */
      anspad_r = xpad_r;
      anspad_i = xpad_i;
    } else {
      ypad_i = helper_i;
      double *yvec = y + j*y_length;
      ypad_r = ans + j*ans_length;
      /* set ypad to zero on the imag part */
      memset(ypad_i, 0, ans_length*sizeof(double));
      /* set ypad to sy at the beginning of the real part */
      memcpy(ypad_r, yvec, y_length*sizeof(double));
      /* set ypad to zero on rest of the real part */
      memset(ypad_r + y_length, 0, (ans_length-y_length)*sizeof(double));
      /* fft of ypad in place */
      fft_factor(ans_length, &maxf, &maxp);
      fft_work(ypad_r, ypad_i, 1, ans_length, 1, 1, work, iwork);
      /* conjugate ypad */
      if (conj) {
        num_elem = ans_length;
        scale_factor = -1;
        increment = 1;
        F77_CALL(dscal)(&num_elem, &scale_factor, ypad_i, &increment);
      }
      anspad_r = ypad_r;
      anspad_i = ypad_i;
    }

    /* There's no blas for complex-complex element wise product */
    /* xpad*ypad: (a + ib) (c + id) = (ac - bd) + i(ad + bc). */
    for (i = 0; i< ans_length; i++) {
      /* use of dtmp  because anspad overwrites either xpad or ypad */
      dtmp = xpad_r[i]*ypad_r[i] - xpad_i[i]*ypad_i[i];
      anspad_i[i] = xpad_r[i]*ypad_i[i] + xpad_i[i]*ypad_r[i];
      anspad_r[i] = dtmp;
    }

    /* inverse of fft of sans */
    fft_factor(ans_length, &maxf, &maxp);
    fft_work(anspad_r, anspad_i, 1, ans_length, 1, -1, work, iwork);
  }

  /* Divide by the length */
  num_elem = ans_length*ans_num_signals;
  scale_factor = 1./ans_length;
  increment = 1;
  F77_CALL(dscal)(&num_elem, &scale_factor, ans, &increment);
  UNPROTECT(3);
  return(sans);
}


/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999-2016   The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/.
 */

 /* Trimmed down and modified by Sergio Oller <sergioller@gmail.com> on
  * Sept 2022
  */

#include <R.h>
#include <Rinternals.h>

#ifndef min
#define min(a, b) ((a < b)?(a):(b))
#define max(a, b) ((a < b)?(b):(a))
#endif

// currently ISNAN includes NAs
#define my_isok(x) (!ISNA(x) & !ISNAN(x))


/* Trimmed down version of cfilter */
SEXP filter(SEXP sx, SEXP sfilter)
{
   if (TYPEOF(sx) != REALSXP || TYPEOF(sfilter) != REALSXP)
       error("invalid input");
    R_xlen_t nx = XLENGTH(sx);
    R_xlen_t nf = XLENGTH(sfilter);

    SEXP ans = allocVector(REALSXP, nx);

    R_xlen_t i, j;
    double z, tmp, *x = REAL(sx), *filter = REAL(sfilter), *out = REAL(ans);

    for (i = 0; i < nf - 1; i++) {
      out[i] = NA_REAL;
    }
    int maxj;
  	for(i = nf; i < nx; i++) {
	    z = 0;
  	  maxj = min(nf, i+1);
	    for(j = max(0, i - nx); j < maxj; j++) {
	      tmp = x[i - j];
	      if(my_isok(tmp)) {
	        z += filter[j] * tmp;
	      } else {
	        out[i] = NA_REAL;
	        goto bad;
	      }
	    }
	    out[i] = z;
    	bad:
	      continue;
    }
    return ans;
}

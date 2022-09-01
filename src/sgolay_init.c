#include "sgolay.h"

static const R_CallMethodDef callMethods[]  = {
  {"cconvolve_circular", (DL_FUNC) &cconvolve_circular, 3},
  {"cfilter", (DL_FUNC) &cfilter, 2},
  {NULL, NULL, 0}
};


void R_init_sgolay(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

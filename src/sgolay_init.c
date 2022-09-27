#include "sgolay.h"

static const R_CallMethodDef callMethods[]  = {
  {"c_pad_and_convolve", (DL_FUNC) &pad_and_convolve, 3},
  {"c_filter", (DL_FUNC) &filter, 2},
  {NULL, NULL, 0}
};


void R_init_sgolay(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

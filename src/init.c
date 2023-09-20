#include "fitLSYS.h"
#include "GRAD_DESC.h"

#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
    {"fitLSYS", (DL_FUNC) fitLSYS, 7},
    {NULL, NULL, 0}
};

void R_init_BGDataExt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

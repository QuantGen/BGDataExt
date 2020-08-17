#define R_NO_REMAP

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callEntries[] = {
    {NULL, NULL, 0}
};

void R_init_BGDataExt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

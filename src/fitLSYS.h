#define R_NO_REMAP

#include <Rinternals.h>

SEXP fitLSYS(SEXP C, SEXP rhs, SEXP b, SEXP active, SEXP RSS, SEXP maxIter, SEXP tolerance);

SEXP fitLSYS(SEXP C, SEXP rhs, SEXP b, SEXP active, SSEXP nIter, SEXP learning_rate);

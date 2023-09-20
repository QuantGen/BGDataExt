#include "GRAD_DESC.h"

SEXP GRAD_DESC(SEXP C, SEXP rhs, SEXP b, SEXP active, SEXP nIter, SEXP learning_rate) {

    int p = Rf_ncols(C);
    R_xlen_t nActive = Rf_xlength(active);
    int nIter = Rf_asInteger(nIter);
    double *pC = REAL(C);
    double *prhs = REAL(rhs);
    b = PROTECT(Rf_duplicate(b));
    double *pb = REAL(b);
    int *pactive = INTEGER(active);


    // Creating a list to return results
    SEXP list = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, b);
    UNPROTECT(2); // b, list
    return list;
}

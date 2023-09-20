#include "GRAD_DESC.h"

SEXP GRAD_DESC(SEXP C, SEXP rhs, SEXP b, SEXP active, SEXP maxIter, SEXP learning_rate) {

    int p = Rf_ncols(C);
    R_xlen_t nActive = Rf_xlength(active);
    int nIter = Rf_asInteger(nIter);
    double *pC = REAL(C);
    double *prhs = REAL(rhs);
    b = PROTECT(Rf_duplicate(b));
    double *pb = REAL(b);
    int *pactive = INTEGER(active);

    for (int iter = 0; iter < nIter; iter++) {
        
        for (int j = 0; j < nActive; j++) { // loop over active predictors
            int k = pactive[j];
            double Ckk = pC[k * (p + 1)];
            double offset = 0;      
            
            for (int m = 0; m < nActive; m++) {
                int n = pactive[m];
                offset += pC[p * k + n] * pb[n];
            }
            
            double gradient=prhs[k]-offset
            double sol=pb[k]-learning_rate*gradient
            pb[k] = sol;
        }
    }

    // Creating a list to return results
    SEXP list = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, b);
    UNPROTECT(2); // b, list
    return list;
}

#include "fitLSYS.h"

SEXP fitLSYS(SEXP C, SEXP rhs, SEXP b, SEXP active, SEXP RSS, SEXP maxIter, SEXP tolerance) {
    int p = Rf_ncols(C);
    R_xlen_t q = Rf_xlength(active);
    int nIter = Rf_asInteger(maxIter);
    double tol = Rf_asReal(tolerance);
    PROTECT(C = Rf_coerceVector(C, REALSXP));
    double *pC = REAL(C);
    PROTECT(rhs = Rf_coerceVector(rhs, REALSXP));
    double *prhs = REAL(rhs);
    PROTECT(b = Rf_coerceVector(b, REALSXP));
    double *pb = REAL(b);
    PROTECT(active = Rf_coerceVector(active, INTSXP));
    int *pactive = INTEGER(active);
    PROTECT(RSS = Rf_coerceVector(RSS, REALSXP));
    double *pRSS = REAL(RSS);
    double RSS0 = pRSS[0] + 0.0;
    int iter = 0;
    while (iter < nIter) {
        iter += 1;
        RSS0 = pRSS[0] + 0.0;
        for (int j = 0; j < q; j++) { // loop over active predictors
            int k = pactive[j];
            double Ckk = pC[k * (p + 1)];
            double offset = 0.0;
            for (int m = 0; m < q; m++) {
                int n = pactive[m];
                offset += pC[p * k + n] * pb[n];
            }
            offset -= Ckk * pb[k];
            double rhs_offset = prhs[k] - offset;
            double sol = rhs_offset / Ckk;
            pRSS[0] += (pow(sol, 2) - pow(pb[k], 2)) * Ckk -2 * (sol - pb[k]) * rhs_offset;
            pb[k] = sol;
        }
        if (((RSS0 - pRSS[0]) / RSS0) < tol) {
            break;
        }
    }
    // Creating a list to return results
    SEXP list = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(list, 0, b);
    SET_VECTOR_ELT(list, 1, RSS);
    UNPROTECT(6);
    return list;
 }

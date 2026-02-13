/*=================================================================
 * dgetrf_mex.c
 * MEX wrapper for LAPACK DGETRF (LU factorization with pivoting)
 *=================================================================*/

#include "mex.h"
#include "lapack.h"

/* Prototype for LAPACK DGETRF */
extern void dgetrf_(const mwSignedIndex *m, const mwSignedIndex *n, 
                    double *A, const mwSignedIndex *lda, mwSignedIndex *ipiv, mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Check inputs */
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("dgetrf_mex:input","One input required (matrix A).");
    }
    if (nlhs > 3) {
        mexErrMsgIdAndTxt("dgetrf_mex:output","Too many output arguments.");
    }

    /* Get input matrix */
    double *A = mxGetPr(prhs[0]);
    mwSignedIndex m = mxGetM(prhs[0]);
    mwSignedIndex n = mxGetN(prhs[0]);
    mwSignedIndex lda = m;

    /* Prepare output matrix (copy of A) */
    plhs[0] = mxDuplicateArray(prhs[0]);
    double *A_out = mxGetPr(plhs[0]);

    /* Prepare pivot array */
    mwSignedIndex min_mn = m < n ? m : n;
    plhs[1] = mxCreateNumericMatrix(min_mn, 1, mxINT32_CLASS, mxREAL);
    mwSignedIndex *ipiv = (mwSignedIndex*) mxGetData(plhs[1]);

    /* Prepare info output */
    mwSignedIndex info = 0;
    if (nlhs > 2) {
        plhs[2] = mxCreateDoubleScalar(0.0); // placeholder, will update later
    }

    /* Call LAPACK DGETRF */
    dgetrf_(&m, &n, A_out, &lda, ipiv, &info);

    /* If requested, return info */
    if (nlhs > 2) {
        mxGetPr(plhs[2])[0] = (double) info;
    }
}

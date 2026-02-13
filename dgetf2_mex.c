#include "mex.h"

/* gfortran / MATLAB LAPACK symbol convention */
extern void dgetf2_(mwSignedIndex *m, mwSignedIndex *n,
                    double *A, mwSignedIndex *lda,
                    mwSignedIndex *ipiv,
                    mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /*
      Usage:
      [A_out, ipiv, info] = dgetf2_mex(A)
    */

    if (nrhs != 1)
        mexErrMsgTxt("Usage: [A_out, ipiv, info] = dgetf2_mex(A)");

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgTxt("A must be a real double matrix.");

    mwSignedIndex m = (mwSignedIndex) mxGetM(prhs[0]);
    mwSignedIndex n = (mwSignedIndex) mxGetN(prhs[0]);
    mwSignedIndex lda = m;
    mwSignedIndex minmn = (m < n) ? m : n;
    mwSignedIndex info = 0;

    /* Duplicate A (factorization is in-place) */
    plhs[0] = mxDuplicateArray(prhs[0]);
    double *A = mxGetPr(plhs[0]);

    /* Create pivot vector */
    plhs[1] = mxCreateNumericMatrix(minmn, 1,
                                    mxINT32_CLASS, mxREAL);
    mwSignedIndex *ipiv =
        (mwSignedIndex *) mxGetData(plhs[1]);

    /* Call LAPACK */
    dgetf2_(&m, &n, A, &lda, ipiv, &info);

    /* Optional INFO output */
    if (nlhs == 3)
        plhs[2] = mxCreateDoubleScalar((double)info);
}

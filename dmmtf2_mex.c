#include "mex.h"

/* gfortran symbol convention */
extern void dmmtf2_(mwSignedIndex *m, mwSignedIndex *n,
                    double *A, mwSignedIndex *lda,
                    double *V, double *W,
                    mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /*
       Usage:
       [A_out, info] = dmmtf2_mex(A, v, w)
    */

    if (nrhs != 3)
        mexErrMsgTxt("Usage: [A_out, info] = dmmtf2_mex(A,v,w)");

    if (nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Get inputs */
    double *v = mxGetPr(prhs[1]);
    double *w = mxGetPr(prhs[2]);

    mwSignedIndex m = (mwSignedIndex) mxGetM(prhs[0]);
    mwSignedIndex n = (mwSignedIndex) mxGetN(prhs[0]);
    mwSignedIndex lda = m;
    mwSignedIndex info = 0;

    /* Duplicate A (factorization is in-place) */
    plhs[0] = mxDuplicateArray(prhs[0]);
    double *A = mxGetPr(plhs[0]);

    /* Call Fortran */
    dmmtf2_(&m, &n, A, &lda,
            v, w,
            &info);

    /* Optional INFO */
    if (nlhs == 2)
        plhs[1] = mxCreateDoubleScalar((double)info);
}

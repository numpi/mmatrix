#include "mex.h"

/* Fortran routine name (lowercase + underscore for gfortran) */
extern void dmmtrff_(mwSignedIndex *m, mwSignedIndex *n,
                     double *A, mwSignedIndex *lda,
                     double *V, double *W,
                     double *WORK,
                     mwSignedIndex *nb,
                     mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Inputs:
       prhs[0] = A
       prhs[1] = v
       prhs[2] = w
       prhs[3] = nb
    */

    if (nrhs != 4)
        mexErrMsgTxt("Usage: [A_out, info] = dmmtrff_mex(A,v,w,nb)");

    if (nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    double *A_in  = mxGetPr(prhs[0]);
    double *v     = mxGetPr(prhs[1]);
    double *w     = mxGetPr(prhs[2]);

    mwSignedIndex nb = (mwSignedIndex) mxGetScalar(prhs[3]);

    mwSignedIndex m = mxGetM(prhs[0]);
    mwSignedIndex n = mxGetN(prhs[0]);
    mwSignedIndex lda = m;
    mwSignedIndex info = 0;

    /* Create output matrix (copy of A) */
    plhs[0] = mxDuplicateArray(prhs[0]);
    double *A = mxGetPr(plhs[0]);

    /* Create work array */
    mwSignedIndex k = (m < n) ? m : n;
    double *work = mxMalloc(k * sizeof(double));

    /* Call Fortran */
    dmmtrff_(&m, &n, A, &lda,
             v, w, work,
             &nb, &info);

    mxFree(work);

    if (nlhs == 2)
        plhs[1] = mxCreateDoubleScalar((double)info);
}

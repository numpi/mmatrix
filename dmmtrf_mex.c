#include "mex.h"

/* Fortran symbol */
extern void dmmtrf_(mwSignedIndex *m, mwSignedIndex *n,
                     double *A, mwSignedIndex *lda,
                     double *V, double *W,
                     double *WORK,
                     mwSignedIndex *nb,
                     mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4)
        mexErrMsgTxt("Usage: [A_out, info] = dmmtrf_mex(A,v,w,nb)");

    if (nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");

    /* Inputs */
    double *A_in = mxGetPr(prhs[0]);
    double *v    = mxGetPr(prhs[1]);
    double *w    = mxGetPr(prhs[2]);
    mwSignedIndex nb = (mwSignedIndex) mxGetScalar(prhs[3]);

    mwSignedIndex m = (mwSignedIndex) mxGetM(prhs[0]);
    mwSignedIndex n = (mwSignedIndex) mxGetN(prhs[0]);
    mwSignedIndex lda = m;
    mwSignedIndex info = 0;

    /* Output matrix */
    plhs[0] = mxDuplicateArray(prhs[0]);
    double *A = mxGetPr(plhs[0]);

    /* Workspace size = n (safer than min(m,n)) */
    double *work = (double*) mxMalloc(n * sizeof(double));

    /* Call Fortran routine */
    dmmtrf_(&m, &n, A, &lda, v, w, work, &nb, &info);

    mxFree(work);

    /* Optional INFO output */
    if (nlhs == 2)
        plhs[1] = mxCreateDoubleScalar((double)info);
}

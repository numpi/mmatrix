#include "mex.h"
#include <string.h>  /* for memcpy */

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

    /* Input sizes */
    mwSignedIndex m  = (mwSignedIndex) mxGetM(prhs[0]);
    mwSignedIndex n  = (mwSignedIndex) mxGetN(prhs[0]);
    mwSignedIndex lda = m;
    mwSignedIndex nb  = (mwSignedIndex) mxGetScalar(prhs[3]);
    mwSignedIndex info = 0;

    /* Duplicate A so original MATLAB matrix is untouched */
    plhs[0] = mxDuplicateArray(prhs[0]);
    double *A = mxGetPr(plhs[0]);

    /* Create PRIVATE COPIES of v and w */
    double *v_in = mxGetPr(prhs[1]);
    double *w_in = mxGetPr(prhs[2]);

    double *v = mxMalloc(n * sizeof(double));
    double *w = mxMalloc(n * sizeof(double));

    memcpy(v, v_in, n * sizeof(double));
    memcpy(w, w_in, n * sizeof(double));

    /* Workspace */
    double *work = mxMalloc(n * sizeof(double));

    /* Call Fortran routine */
    dmmtrf_(&m, &n, A, &lda,
            v, w, work,
            &nb, &info);

    /* Free allocated memory */
    mxFree(v);
    mxFree(w);
    mxFree(work);

    if (nlhs == 2)
        plhs[1] = mxCreateDoubleScalar((double)info);
}

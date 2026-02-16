#include "mex.h"
#include <string.h>

/* gfortran symbol convention */
extern void dmmtf2_(mwSignedIndex *m, mwSignedIndex *n,
                    double *A, mwSignedIndex *lda,
                    double *V, double *W,
                    mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* --------------------------------------------------
       Argument checking
    -------------------------------------------------- */
    if (nrhs != 3)
        mexErrMsgIdAndTxt("dmmtf2:nrhs",
                          "Usage: [A_out, info] = dmmtf2_mex(A,v,w)");

    if (nlhs > 2)
        mexErrMsgIdAndTxt("dmmtf2:nlhs",
                          "Too many output arguments.");

    if (!mxIsDouble(prhs[0]) || mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt("dmmtf2:A",
                          "A must be full double matrix.");

    if (!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]))
        mexErrMsgIdAndTxt("dmmtf2:vw",
                          "v and w must be double vectors.");

    /* --------------------------------------------------
       Dimensions
    -------------------------------------------------- */
    mwSize m = mxGetM(prhs[0]);
    mwSize n = mxGetN(prhs[0]);

    if (m != n)
        mexErrMsgIdAndTxt("dmmtf2:square",
                          "A must be square.");

    mwSize vn = mxGetNumberOfElements(prhs[1]);
    mwSize wn = mxGetNumberOfElements(prhs[2]);

    if (vn != n || wn != n)
        mexErrMsgIdAndTxt("dmmtf2:dimension",
                          "v and w must have length n.");

    mwSignedIndex M = (mwSignedIndex)m;
    mwSignedIndex N = (mwSignedIndex)n;
    mwSignedIndex lda = M;
    mwSignedIndex info = 0;

    /* --------------------------------------------------
       Duplicate A (factorization is in-place)
    -------------------------------------------------- */
    plhs[0] = mxDuplicateArray(prhs[0]);
    double *A = mxGetPr(plhs[0]);

    /* --------------------------------------------------
       Protect v and w (copy them)
    -------------------------------------------------- */
    const double *v_in = mxGetPr(prhs[1]);
    const double *w_in = mxGetPr(prhs[2]);

    double *v = mxMalloc(n * sizeof(double));
    double *w = mxMalloc(n * sizeof(double));

    memcpy(v, v_in, n * sizeof(double));
    memcpy(w, w_in, n * sizeof(double));

    /* --------------------------------------------------
       Call Fortran
    -------------------------------------------------- */
    dmmtf2_(&M, &N, A, &lda,
             v, w,
            &info);

    /* --------------------------------------------------
       Free memory
    -------------------------------------------------- */
    mxFree(v);
    mxFree(w);

    if (nlhs == 2)
        plhs[1] = mxCreateDoubleScalar((double)info);
}

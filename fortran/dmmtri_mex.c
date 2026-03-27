#include "mex.h"

/* Fortran routine (blocked M-matrix inverse) */
extern void dmmtri_(mwSignedIndex *n,
                    double *A,
                    mwSignedIndex *lda,
                    double *work,
                    mwSignedIndex *lwork,
                    mwSignedIndex *info);

/* MEX gateway for dmmtri */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1)
        mexErrMsgTxt("Usage: Ainv = dmmtri_mex(A)");
    if (nlhs > 2)
        mexErrMsgTxt("Too many output arguments");

    double *A_in = mxGetPr(prhs[0]);
    mwSignedIndex n = (mwSignedIndex) mxGetM(prhs[0]);
    mwSignedIndex lda = n;
    mwSignedIndex info;

    /* Copy input matrix to output */
    plhs[0] = mxDuplicateArray(prhs[0]);
    double *A = mxGetPr(plhs[0]);

    /* --- Workspace query --- */
    mwSignedIndex lwork = -1;
    double work_query;
    dmmtri_(&n, A, &lda, &work_query, &lwork, &info);

    if (info != 0) {
        mexErrMsgTxt("dmmtri workspace query failed: invalid input.");
    }

    /* Allocate optimal workspace */
    lwork = (mwSignedIndex) work_query;
    double *work = mxMalloc(lwork * sizeof(double));

    /* Call Fortran routine with optimal workspace */
    dmmtri_(&n, A, &lda, work, &lwork, &info);

    /* Free workspace */
    mxFree(work);

    /* Optional second output: info */
    if (nlhs == 2)
        plhs[1] = mxCreateDoubleScalar((double)info);

    if (info != 0)
        mexWarnMsgTxt("dmmtri: matrix may be singular or ill-conditioned");
}
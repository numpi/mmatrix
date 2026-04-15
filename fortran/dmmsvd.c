#include "mex.h"
#include <math.h>

/* Fortran routines */
extern void dmmltf_(mwSignedIndex *n,
                     double *A,
                     mwSignedIndex *lda,
                     double *s,
                     mwSignedIndex *piv,
                     mwSignedIndex *info);

extern void dmmsvd_(mwSignedIndex *n,
                     double *L, mwSignedIndex *ldl,
                     double *D,
                     double *U, mwSignedIndex *ldu,
                     double *S,
                     double *UMAT, mwSignedIndex *lduf,
                     double *VMAT, mwSignedIndex *ldvf,
                     double *work, mwSignedIndex *lwork,
                     mwSignedIndex *iwork,
                     mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* ================= Input check ================= */
    if (nrhs != 2)
        mexErrMsgTxt("Usage: [U,S,V] = dmmsvd(Aoff, s)");

    mwSignedIndex n = mxGetM(prhs[0]);

    if (mxGetN(prhs[0]) != n)
        mexErrMsgTxt("Aoff must be square");

    if (mxGetM(prhs[1]) != n)
        mexErrMsgTxt("s must be a vector of length n");

    /* Duplicate inputs (since DMMLTF overwrites) */
    mxArray *Ain = mxDuplicateArray(prhs[0]);
    mxArray *sin = mxDuplicateArray(prhs[1]);

    double *A = mxGetPr(Ain);
    double *s = mxGetPr(sin);

    mwSignedIndex *piv =
        (mwSignedIndex*)mxMalloc(n*sizeof(mwSignedIndex));

    mwSignedIndex info;

    /* ================= DMMLTF ================= */
    dmmltf_(&n, A, &n, s, piv, &info);

    if (info != 0)
        mexErrMsgTxt("DMMLTF failed");

    /* ================= Extract L, D, U ================= */
    double *L = (double*)mxMalloc(n*n*sizeof(double));
    double *U = (double*)mxMalloc(n*n*sizeof(double));
    double *D = (double*)mxMalloc(n*sizeof(double));

    for (mwSignedIndex i = 0; i < n; i++) {
        D[i] = A[i + i*n];

        for (mwSignedIndex j = 0; j < n; j++) {
            if (i > j) {
                L[i + j*n] = A[i + j*n];
                U[i + j*n] = 0.0;
            }
            else if (i == j) {
                L[i + j*n] = 1.0;
                U[i + j*n] = 1.0;
            }
            else {
                L[i + j*n] = 0.0;
                U[i + j*n] = A[i + j*n] / A[i + i*n];
            }
        }
    }

    /* ================= Outputs ================= */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL); /* U */
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL); /* S */
    plhs[2] = mxCreateDoubleMatrix(n, n, mxREAL); /* V */

    double *Uout = mxGetPr(plhs[0]);
    double *Sout = mxGetPr(plhs[1]);
    double *Vout = mxGetPr(plhs[2]);

    /* ================= Workspace ================= */
    mwSignedIndex lwork = n*n + 10*n;

    double *work = (double*)mxMalloc(lwork*sizeof(double));

    mwSignedIndex *iwork =  (mwSignedIndex*)mxMalloc((n+10)*sizeof(mwSignedIndex));

    /* ================= DMMSVD ================= */
    dmmsvd_(&n, L, &n, D, U, &n,
             Sout,
             Uout, &n,
             Vout, &n,
             work, &lwork,
             iwork,
             &info);

    if (info != 0)
        mexErrMsgTxt("DMMSVD failed");

    /* ================= Cleanup ================= */
    mxFree(L);
    mxFree(U);
    mxFree(D);
    mxFree(work);
    mxFree(iwork);
    mxFree(piv);

    mxDestroyArray(Ain);
    mxDestroyArray(sin);
}
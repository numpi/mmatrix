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
                     mwSignedIndex *piv,    
                     mwSignedIndex *info);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
        mexErrMsgTxt("Usage: [U,S,V] = dmmsvd_mex(Aoff, s)");

    mwSignedIndex n = mxGetM(prhs[0]);

    if (mxGetN(prhs[0]) != n)
        mexErrMsgTxt("Aoff must be square");

    /* Duplicate inputs (Fortran will overwrite) */
    mxArray *Ain = mxDuplicateArray(prhs[0]);
    mxArray *sin = mxDuplicateArray(prhs[1]);

    double *A = mxGetPr(Ain);
    double *s = mxGetPr(sin);

    /* pivot array */
    mwSignedIndex *piv =
        (mwSignedIndex*)mxMalloc(n * sizeof(mwSignedIndex));

    mwSignedIndex info;

/* ============================================================
   DMMLTF
   ============================================================ */
    dmmltf_(&n, A, &n, s, piv, &info);

    if (info != 0)
        mexErrMsgTxt("DMMLTF failed");

/* ============================================================
   Extract L, D, U from A
   ============================================================ */

    double *L = (double*)mxCalloc(n * n, sizeof(double));
    double *U = (double*)mxCalloc(n * n, sizeof(double));
    double *D = (double*)mxCalloc(n, sizeof(double));

    for (mwSignedIndex i = 0; i < n; i++)
    {
        D[i] = A[i + i * n];

        for (mwSignedIndex j = 0; j < n; j++)
        {
            if (i > j)
            {
                L[i + j * n] = A[i + j * n];
                U[i + j * n] = 0.0;
            }
            else if (i == j)
            {
                L[i + j * n] = 1.0;
                U[i + j * n] = 1.0;
            }
            else
            {
                L[i + j * n] = 0.0;
                U[i + j * n] = A[i + j * n] / A[i + i * n];
            }
        }
    }

/* ============================================================
   Outputs
   ============================================================ */

    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n, n, mxREAL);

    double *Uout = mxGetPr(plhs[0]);
    double *Sout = mxGetPr(plhs[1]);
    double *Vout = mxGetPr(plhs[2]);

/* ============================================================
    Workspace
   ============================================================ */

    mwSignedIndex lwork = n*n + 10*n;

    double *work =
        (double*)mxMalloc(lwork * sizeof(double));

    mwSignedIndex *iwork =
        (mwSignedIndex*)mxMalloc((n + 10) * sizeof(mwSignedIndex));

/* ============================================================
    DMMSVD  
   ============================================================ */

    dmmsvd_(&n,
            L, &n,
            D,
            U, &n,
            Sout,
            Uout, &n,
            Vout, &n,
            work, &lwork,
            iwork,
            piv,   
            &info);

    if (info != 0)
        mexErrMsgTxt("DMMSVD failed");

/* ============================================================
   cleanup
   ============================================================ */

    mxFree(L);
    mxFree(U);
    mxFree(D);
    mxFree(work);
    mxFree(iwork);
    mxFree(piv);

    mxDestroyArray(Ain);
    mxDestroyArray(sin);
}

/* --------------------------------------------------
   Generic Lapack factorization Mex function wrapper

   This code contains a Mex wrapper around a Fortran subroutine
   specified in the macro FUNNAME

   If SOLVE is defined, the wrapper also solves a linear system.
   If LEFT_SOLVE, this linear system is X A = B instead of A X = B.
   
   If HAS_WORK is provided, a temporary array WORK is allocated
   and passed to the subroutine. 
   This array normally has size N, but if WORK_SIZE_NB, then it has
   size NB, and if WORK_SIZE_DMMTRFF then the size is either 2*NB or N,
   according to what DMMTRFF needs.

   Inputs:

   Mandatory:
   1 (A)
   +2 (v,w) if HAS_TRIPLET
   +1 (b) if SOLVE
   Optional:
   +1 (NB) if HAS_NB

   Outputs:
   1 (x) if SOLVE
   +1 (A)
   +1 (ipiv) if HAS_IPIV
   +1 (info)
   +1 (w) if HAS_TRIPLET

-------------------------------------------------- */


#include "mex.h"

#ifndef FUNNAME
#error "FUNNAME must be defined at compile time, e.g. -DFUNNAME=dmmtf2"
#endif

/* Macro helpers */
#define PASTE2(a,b) a##b
#define EXPAND_AND_PASTE2(a,b) PASTE2(a,b)
#define F77_NAME(name) EXPAND_AND_PASTE2(name, _)

#define STR1(x) #x
#define STR(x)  STR1(x)
#define FUNNAME_STR STR(FUNNAME)

/* gfortran symbol convention: append underscore */
extern void F77_NAME(FUNNAME)(mwSignedIndex *m, mwSignedIndex *n,
                              double *A, mwSignedIndex *lda,
#ifdef HAS_TRIPLET                              
                              double *V, double *W,
#endif                              
#ifdef HAS_WORK                              
                              double *WORK,
#endif                             
#ifdef HAS_NB
                              mwSignedIndex *nb,
#endif       
#ifdef HAS_IPIV
                              mwSignedIndex *ipiv,
#endif                              
                              mwSignedIndex *info);

/* LAPACK function signatures */
extern mwSignedIndex F77_NAME(ilaenv)(mwSignedIndex *ispec,
                              const char *name,
                              const char *opts,
                              mwSignedIndex *n1,
                              mwSignedIndex *n2,
                              mwSignedIndex *n3,
                              mwSignedIndex *n4);

extern void F77_NAME(dlaswp)(mwSignedIndex *n, double *A, mwSignedIndex *lda,
                             mwSignedIndex *k1, mwSignedIndex *k2,
                             mwSignedIndex *ipiv, mwSignedIndex *incx);

extern void F77_NAME(dtrsm)(const char *side, const char *uplo, const char *transa,
                             const char *diag, mwSignedIndex *m, mwSignedIndex *n,
                             const double *alpha, const double *A, mwSignedIndex *lda,
                             double *B, mwSignedIndex *ldb);


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* --------------------------------------------------
       Argument checking
    -------------------------------------------------- */
    int min_inputs = 1;
    int optional_inputs = 0;
    int b_position = 1; /*position of input b (if it exists); 1 = second input*/
#ifdef HAS_TRIPLET
    min_inputs += 2;
    b_position += 2;
#endif
#ifdef SOLVE
    min_inputs += 1;
#endif
#ifdef HAS_NB
    optional_inputs += 1;
#endif

    int max_outputs = 2; /*A and info*/ 
#ifdef SOLVE
    max_outputs += 1;
#endif
#ifdef HAS_IPIV
    max_outputs += 1;
#endif
#ifdef HAS_TRIPLET
    max_outputs += 1;
#endif
 
    if ((nrhs < min_inputs ) || (nrhs > min_inputs + optional_inputs))
        mexErrMsgIdAndTxt(FUNNAME_STR ":nrhs",
                          "Wrong number of inputs");

    if (nlhs > max_outputs)
        mexErrMsgIdAndTxt(FUNNAME_STR ":nlhs",
                          "Too many output arguments.");

    if (!mxIsDouble(prhs[0]) || mxIsSparse(prhs[0]))
        mexErrMsgIdAndTxt(FUNNAME_STR ":A",
                          "A must be a full double matrix.");
#ifdef HAS_TRIPLET
    if (!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]))
        mexErrMsgIdAndTxt(FUNNAME_STR ":vw",
                          "v and w must be double vectors.");
#endif
#ifdef SOLVE
    if (!mxIsDouble(prhs[b_position]))
        mexErrMsgIdAndTxt(FUNNAME_STR ":b",
                          "b must be a double vector or matrix.");
#endif
    if (mxIsComplex(prhs[0]) 
#ifdef HAS_TRIPLET        
        || mxIsComplex(prhs[1]) || mxIsComplex(prhs[2])
#endif
#ifdef SOLVE
        || mxIsComplex(prhs[b_position])
#endif        
       )
        mexErrMsgIdAndTxt(FUNNAME_STR ":complex",
                          "Inputs must be real (non-complex).");

    /* --------------------------------------------------
       Dimensions
    -------------------------------------------------- */
    mwSignedIndex one = 1;
    mwSignedIndex minus_one = -1;
    mwSize m = mxGetM(prhs[0]);
    mwSize n = mxGetN(prhs[0]);
    mwSize minmn = m < n ? m : n;
    mwSignedIndex M   = (mwSignedIndex)m; /*Lapack wants these with a different type (signed int)*/
    mwSignedIndex N   = (mwSignedIndex)n; /*Uppercase = signed version*/
    mwSignedIndex LDA = M;

#ifdef HAS_TRIPLET
#ifdef LEFT_TRIPLET
    mwSize vlength = m;
    mwSize wlength = n;
#else
    mwSize vlength = n;
    mwSize wlength = m;
#endif   /*LEFT_TRIPLET*/
    if (mxGetNumberOfElements(prhs[1]) != vlength)
        mexErrMsgIdAndTxt(FUNNAME_STR ":dimension",
                          "v has the wrong length.");
    if (mxGetNumberOfElements(prhs[2]) != wlength)
        mexErrMsgIdAndTxt(FUNNAME_STR ":dimension",
                          "w has the wrong length.");
#endif  /*HAS_TRIPLET*/

#ifdef SOLVE
#ifdef LEFT_SOLVE
    mwSize sizerhs = mxGetM(prhs[b_position]);
    if (mxGetN(prhs[b_position]) != m)
        mexErrMsgIdAndTxt(FUNNAME_STR ":dimension",
                          "wrong size(b, 2).");
    mwSignedIndex SIZERHS   = (mwSignedIndex)sizerhs;
    mwSignedIndex LDB = SIZERHS;
#else /*not LEFT_SOLVE*/
    mwSize sizerhs = mxGetN(prhs[b_position]);
    if (mxGetM(prhs[b_position]) != n)
        mexErrMsgIdAndTxt(FUNNAME_STR ":dimension",
                          "wrong size(b, 1).");
    mwSignedIndex SIZERHS   = (mwSignedIndex)sizerhs;
    mwSignedIndex LDB = N;
#endif /*LEFT_SOLVE*/
#endif /*SOLVE*/

    mwSignedIndex worksize = (mwSignedIndex)minmn;
    
#ifdef HAS_NB
    mwSignedIndex nb;
    
    if (nrhs > min_inputs)
        nb  = (mwSignedIndex) mxGetScalar(prhs[min_inputs]);
    else
        nb = F77_NAME(ilaenv)(&one, "DGETRF", " ", &M, &N, &minus_one, &minus_one);

    worksize = nb;
#endif
#ifdef WORK_SIZE_DMMTRFF
    if(nb == 1)
        worksize = n;
    else
        worksize = 2*nb;
#endif
    mwSignedIndex info = 0;

    /* --------------------------------------------------
       Duplicate A (factorization is in-place)
    -------------------------------------------------- */
    mxArray *Ap = mxDuplicateArray(prhs[0]);
    double *A = mxGetPr(Ap);

#ifdef HAS_TRIPLET
    /* --------------------------------------------------
       v is read-only: no copy/duplication
    -------------------------------------------------- */
    const double *v = mxGetPr(prhs[1]);
    /* --------------------------------------------------
       Duplicate w (FUNNAME_ overwrites W)
    -------------------------------------------------- */
    mxArray *wp = mxDuplicateArray(prhs[2]);
    double *w = mxGetPr(wp);    
#endif
#ifdef SOLVE
    /* --------------------------------------------------
       Duplicate b (we will overwrites b with x)
    -------------------------------------------------- */
    mxArray *bp = mxDuplicateArray(prhs[b_position]);
    double *b = mxGetPr(bp);
#endif

#ifdef HAS_WORK    
    /* Workspace */
    double *work = mxMalloc(worksize * sizeof(double));
#endif
    
#ifdef HAS_IPIV
    /* Prepare pivot array */
    mwSignedIndex min_mn = m < n ? m : n;
    mxArray *ipivp = mxCreateNumericMatrix(min_mn, 1, mxINT64_CLASS, mxREAL);
    mwSignedIndex *ipiv = (mwSignedIndex*) mxGetData(ipivp);
#endif

    /* --------------------------------------------------
       Call Fortran
    -------------------------------------------------- */
    F77_NAME(FUNNAME)(&M, &N, A, &LDA,
#ifdef HAS_TRIPLET                      
                      (double *)v, w,
#endif                      
#ifdef HAS_WORK                              
                              work,
#endif                             
#ifdef HAS_NB
                              &nb,
#endif       
#ifdef HAS_IPIV
                              ipiv,
#endif 
                      &info);

#ifdef SOLVE
    double oneD0 = 1.0;
    /* --------------------------------------------------
       Use factorization to solve the linear system
    -------------------------------------------------- */
#ifdef LEFT_SOLVE
#ifdef HAS_IPIV
    #error "LEFT_SOLVE with HAS_IPIV is not implemented"
#endif /*HAS_IPIV*/
    F77_NAME(dtrsm)("R", "U", "N", "N",
                    &SIZERHS, &M, &oneD0, A, &LDA, b, &LDB);
    
    F77_NAME(dtrsm)("R", "L", "N", "U",
                    &SIZERHS, &M, &oneD0, A, &LDA, b, &LDB);
#else /*not LEFT_SOLVE*/
#ifdef HAS_IPIV
    F77_NAME(dlaswp)(&SIZERHS, b, &LDB, &one, &minmn, ipiv, &one);
#endif /*HAS_IPIV*/
    F77_NAME(dtrsm)("L", "L", "N", "U",
                    &N, &SIZERHS, &oneD0, A, &LDA, b, &LDB);

    F77_NAME(dtrsm)("L", "U", "N", "N",
                    &N, &SIZERHS, &oneD0, A, &LDA, b, &LDB);
#endif /*LEFT_SOLVE*/
#endif /*SOLVE*/

    /* --------------------------------------------------
       Cleanup used arrays, apart from the ones
       that are returned.
    -------------------------------------------------- */

#ifdef HAS_WORK
    if (work) { mxFree(work); work = NULL; }
#endif

    int filled_outputs = 0;
 
#ifdef SOLVE
    if (nlhs > filled_outputs) {
        plhs[filled_outputs] = bp;
        filled_outputs += 1;
    }
    else
        mxDestroyArray(bp);
#endif

    if (nlhs > filled_outputs) {
        plhs[filled_outputs] = Ap;
        filled_outputs += 1;
    }
    else
        mxDestroyArray(Ap);

#ifdef HAS_IPIV
    if (nlhs > filled_outputs){
        plhs[filled_outputs] = ipivp;
        filled_outputs += 1;
    }
    else
        mxDestroyArray(ipivp);
#endif

    if (nlhs > filled_outputs){
        plhs[filled_outputs] = mxCreateDoubleScalar((double)info);
        filled_outputs += 1;
    }

#ifdef HAS_TRIPLET
    if (nlhs >= filled_outputs){
        plhs[filled_outputs] = wp;
        filled_outputs += 1;
    }
    else
        mxDestroyArray(wp);
#endif    
}

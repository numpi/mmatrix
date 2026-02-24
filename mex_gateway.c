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

/* --------------------------------------------------
   Generic Lapack factorization function wrapper
   
   inputs:
   1 numeric (A) or 3 (A,v,w) if HAS_TRIPLET
   followed by non-numeric NB if HAS_NB

   outputs:
   1 numeric (A) or 2 (A,ipiv) if HAS_IPIV
   followed by info (optional) and w (optional) if HAS_TRIPLET

-------------------------------------------------- */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* --------------------------------------------------
       Argument checking
    -------------------------------------------------- */
    int numinputs = 1; /*numeric inputs: A or A,v,w */
#ifdef HAS_TRIPLET
    numinputs += 2;
#endif
    int nonnuminputs = 0; /*non-numeric inputs: nb*/
#ifdef HAS_NB
    nonnuminputs += 1;
#endif

    int numoutputs = 1; /*number of numeric outputs (A and maybe ipiv) */
#ifdef HAS_IPIV
    numoutputs += 1;
#endif
    int nonnumoutputs = 1; /*non-numeric outputs: info, w*/
#ifdef HAS_TRIPLET
    nonnumoutputs += 1;
#endif
 
    if (nrhs != numinputs+nonnuminputs)
        mexErrMsgIdAndTxt(FUNNAME_STR ":nrhs",
                          "Wrong number of inputs");

    if (nlhs > numoutputs+nonnumoutputs)
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
    if (mxIsComplex(prhs[0]) 
#ifdef HAS_TRIPLET        
        || mxIsComplex(prhs[1]) || mxIsComplex(prhs[2])
#endif        
       )
        mexErrMsgIdAndTxt(FUNNAME_STR ":complex",
                          "Inputs must be real (non-complex).");

    /* --------------------------------------------------
       Dimensions
    -------------------------------------------------- */
    mwSize m = mxGetM(prhs[0]);
    mwSize n = mxGetN(prhs[0]);
    mwSize minmn = m < n ? m : n;

    if (m != n)
        mexErrMsgIdAndTxt(FUNNAME_STR ":square",
                          "A must be square.");

#ifdef HAS_TRIPLET    
    if (mxGetNumberOfElements(prhs[1]) != n ||
        mxGetNumberOfElements(prhs[2]) != n)
        mexErrMsgIdAndTxt(FUNNAME_STR ":dimension",
                          "v and w must have length n.");
#endif    

    mwSignedIndex M   = (mwSignedIndex)m;
    mwSignedIndex N   = (mwSignedIndex)n;
    mwSignedIndex worksize = (mwSignedIndex)minmn;
    mwSignedIndex lda = M;
#ifdef HAS_NB
    mwSignedIndex nb  = (mwSignedIndex) mxGetScalar(prhs[3]);
    worksize = nb;
#endif
    mwSignedIndex info = 0;

    /* --------------------------------------------------
       Duplicate A (factorization is in-place)
    -------------------------------------------------- */
    plhs[0] = mxDuplicateArray(prhs[0]);
    double *A = mxGetPr(plhs[0]);

#ifdef HAS_TRIPLET    
    /* --------------------------------------------------
       v is read-only: no copy/duplication
    -------------------------------------------------- */
    const double *v = mxGetPr(prhs[1]);
    /* --------------------------------------------------
       Duplicate w (FUNNAME_ overwrites W)
    -------------------------------------------------- */
    plhs[numoutputs+1] = mxDuplicateArray(prhs[2]);
    double *w = mxGetPr(plhs[numoutputs+1]);    
#endif    

#ifdef HAS_WORK    
    /* Workspace */
    double *work = mxMalloc(worksize * sizeof(double));
#endif
    
#ifdef HAS_IPIV
    /* Prepare pivot array */
    mwSignedIndex min_mn = m < n ? m : n;
    plhs[1] = mxCreateNumericMatrix(min_mn, 1, mxINT64_CLASS, mxREAL);
    mwSignedIndex *ipiv = (mwSignedIndex*) mxGetData(plhs[1]);
#endif

    /* --------------------------------------------------
       Call Fortran
    -------------------------------------------------- */
    F77_NAME(FUNNAME)(&M, &N, A, &lda,
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


#ifdef HAS_WORK
    /* --------------------------------------------------
       Cleanup work array
    -------------------------------------------------- */
    if (work) { mxFree(work); work = NULL; }
#endif

    /* --------------------------------------------------
       info output if requested
    -------------------------------------------------- */
    if (nlhs >= numoutputs+1)
        plhs[numoutputs] = mxCreateDoubleScalar((double)info);

#ifdef HAS_TRIPLET   
    /* --------------------------------------------------
       If w is not requested as output, free it and clear plhs[2]
    -------------------------------------------------- */
    if (nlhs < numoutputs + 2) {
        mxDestroyArray(plhs[numoutputs+1]);
        plhs[numoutputs+1] = NULL;
    }
#endif    
}

"""
This script tests:
1. GTH algorithm implemented in Fortran (via DMMTF2 and DMMTRF).
2. Standard LAPACK routines (DGETF2 and DGETRF) for LU factorization.
3. A pure Python implementation of GTH.
"""
import numpy as np
# np.__config__.show()

import time
import mpmath as mp
import scipy.linalg as la
from scipy.linalg.lapack import get_lapack_funcs
import os
import sys
import subprocess
# ------------------------------------------------------------
# Compile Fortran modules (once)
def compile_fortran_modules():
    # First, check if dmmtf2 exists as a separate module
    #if not os.path.exists("dmmtf2.so"):
        cmd = [sys.executable, "-m", "numpy.f2py", "-c", "dmmtf2.f",
               "-m", "dmmtf2", "--fcompiler=gnu95", "-lopenblas"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            print("dmmtf2 compilation failed:")
            print(res.stderr)
    # if not os.path.exists("dmmtrf2.so"):
    #     cmd = [sys.executable, "-m", "numpy.f2py", "-c", "dmmtrf2.f",
    #            "-m", "dmmtrf2", "--fcompiler=gnu95", "-lopenblas"]
    #     res = subprocess.run(cmd, capture_output=True, text=True)
    #     if res.returncode != 0:
    #         print("dmmtrf2 compilation failed:")
    #         print(res.stderr)
    # Compile dgetf2
    #if not os.path.exists("dgetf2.so"):
        cmd = [sys.executable, "-m", "numpy.f2py", "-c", "dgetf2.f",
               "-m", "dgetf2", "--fcompiler=gnu95", "-lopenblas"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            print("dgetf2 compilation failed:")
            print(res.stderr)
            sys.exit(1)
    
    # Compile dmmtrf WITH dmmtf2 included
    #if not os.path.exists("dmmtrf.so"):
        # Method 1: Compile both files together
        cmd = [sys.executable, "-m", "numpy.f2py", "-c", 
               "dmmtf2.f", "dmmtrf.f",  # Include both source files
               "-m", "dmmtrf", 
               "--fcompiler=gnu95", 
               "-lopenblas"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            print("dmmtrf compilation failed (method 1):")
            print(res.stderr)
            
            # Method 2: Try with manual gfortran compilation
            print("\nTrying manual gfortran compilation...")
            cmd = ["gfortran", "-fPIC", "-shared", 
                   "-o", "dmmtrf.so", 
                   "dmmtf2.f", "dmmtrf.f",  # Both files
                   "-lopenblas", "-llapack"]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                print("Manual compilation failed:")
                print(res.stderr)
                sys.exit(1)
# ------------------------------------------------------------
# GTH Fortran  
def gth_fortran(A, v, w, module):
    m, n = A.shape
    A_f = np.asfortranarray(A)
    info = np.array(0, dtype=np.int32)
    module.dmmtf2(m=m, n=n, a=A_f, lda=m, v=v, w=w, info=info)
    return A_f

# ------------------------------------------------------------
# DGETF2 Fortran  
def dgetf2_fortran(A, module):
    m, n = A.shape
    A_f = np.asfortranarray(A.copy())
    ipiv = np.zeros(min(m, n), dtype=np.int32)
    info = np.array(0, dtype=np.int32)
    module.dgetf2(m=m, n=n, a=A_f, lda=m, ipiv=ipiv, info=info)
    return A_f, ipiv, info[()]
# ------------------------------------------------------------
# DMMTRF Fortran 
def dmmtrf_fortran(A, v, w, module, nb):
    m, n = A.shape
    A_f = np.asfortranarray(A)
    info = np.array(0, dtype=np.int32)
    nb_val = np.array(nb, dtype=np.int32)
    work = np.zeros(nb, dtype=np.float64, order="F")

    module.dmmtrf(m=m, n=n, a=A_f, lda=m, v=v,  w=w, work=work, nb=nb_val, info=info)
    return A_f

# ------------------------------------------------------------
# GTH Python  
def gth_python(A, v, w):
    n = A.shape[0]
    L = np.eye(n)
    U = A.copy()
    for k in range(n - 1):
        pivot = (w[k] - U[k, k+1:] @ v[k+1:]) / v[k]
        U[k, k] = pivot
        L[k+1:, k] = U[k+1:, k] / pivot
        U[k+1:, k] = 0.0
        w[k+1:] -= L[k+1:, k] * w[k]
        U[k+1:, k+1:] -= np.outer(L[k+1:, k], U[k, k+1:])
    U[-1, -1] = w[-1] / v[-1]
    return L, U
# ------------------------------------------------------------
def m_matrix(n, seed):
    rng = np.random.default_rng(seed)
    A = -rng.uniform(0, 1e-1, (n, n))
    np.fill_diagonal(A, 0)
    diag_values = np.sum(np.abs(A), axis=1) + 1.0
    np.fill_diagonal(A, diag_values)
    return A
# ------------------------------------------------------------
# compute w = A v exactly
def exact_rhs(A):
    n = A.shape[0]
    mp.mp.dps = 100  # high precision
    v_mp = mp.ones(n, 1)
    A_mp = mp.matrix(A.tolist())
    w_mp = A_mp * v_mp
    w = np.array(w_mp.tolist(), dtype=np.float64).ravel()
    v = np.ones(n, dtype=np.float64)
    #v = np.random.rand(n)
    return v, w
# ------------------------------------------------------------
# Solve LUx = b
def solve_lu(L, U, b):
    y = la.solve_triangular(L, b, lower=True, check_finite=False)
    return la.solve_triangular(U, y, check_finite=False)
# ------------------------------------------------------------
# Relative error
def rel_error(x, x_exact):
    return np.max(np.abs(x - x_exact) / np.abs(x_exact))
    #return la.norm(x - x_exact, np.inf) / la.norm(x_exact, np.inf)
# ------------------------------------------------------------
# Main test
def run_test(n_list, num_repeats):
    import dmmtf2
    import dgetf2
    import dmmtrf
#    import dmmtrf2

    print(f"{'n':>6} "
          f"{'dmmtf2 time':>12} {'dmmtf2 err':>12} "
          f"{'dgetf2 time':>12} {'dgetf2 err':>12} "
          f"{'dmmtrf time':>12} {'dmmtrf err':>12} "
          f"{'dgetrf time':>12} {'dgetrf err':>12} "
          f"{'GTH-P time':>12} {'GTH-P err':>12}")
    print("-" * 160)

    for n in n_list:
        A = m_matrix(n, 1)
        nb = 128  # block size for dmmtrf
        v, w = exact_rhs(A)
        b = w.copy()
        x_exact = np.ones(n)

        # Average times over multiple runs
        t_gth_f, t_dgetf2, t_dmmtrf, t_lapack, t_gth_p = 0, 0, 0, 0, 0
        err_gth_f, err_dgetf2, err_dmmtrf, err_lapack, err_gth_p = 0, 0, 0, 0, 0
        
        for _ in range(num_repeats):
            # Fortran GTH (dmmtf2)
            t0 = time.perf_counter()
            A_f = gth_fortran(A, v.copy(), w.copy(), dmmtf2)
            t_gth_f += (time.perf_counter() - t0)
            Lf = np.tril(A_f, -1) + np.eye(n)
            Uf = np.triu(A_f)
            x_gth_f = solve_lu(Lf, Uf, b)
            err_gth_f = rel_error(x_gth_f, x_exact)
            
            # Fortran DGETF2
            t0 = time.perf_counter()
            lu_f2, ipiv_f2, _ = dgetf2_fortran(A, dgetf2)
            t_dgetf2 += (time.perf_counter() - t0) 
            L = np.tril(lu_f2, -1) + np.eye(n)
            U = np.triu(lu_f2)
            b_perm = b.copy()
            for i in range(len(ipiv_f2)):
                if ipiv_f2[i] != i + 1:  
                    b_perm[[i, ipiv_f2[i]-1]] = b_perm[[ipiv_f2[i]-1, i]]
            y = la.solve_triangular(L, b_perm, lower=True, check_finite=False)
            x_dgetf2 = la.solve_triangular(U, y, check_finite=False)
            err_dgetf2 = rel_error(x_dgetf2, x_exact)
            
            # Fortran DMMTRF (same interface as dmmtf2)
            t0 = time.perf_counter()
            A_trfa = dmmtrf_fortran(A, v.copy(), w.copy(), dmmtrf, nb)
            t_dmmtrf += (time.perf_counter() - t0)
            L = np.tril(A_trfa, -1) + np.eye(n)
            U = np.triu(A_trfa)
            x_dmmtrf = solve_lu(L, U, b)
            err_dmmtrf = rel_error(x_dmmtrf, x_exact)

            # LAPACK DGETRF
            A_f = np.asfortranarray(A)
            dgetrf = get_lapack_funcs('getrf', (A_f,))
            t0 = time.perf_counter()
            lu_rf, piv_rf, _ = dgetrf(A_f)
            t_lapack += (time.perf_counter() - t0)
            x_lapack = la.lu_solve((lu_rf, piv_rf), b)
            err_lapack = rel_error(x_lapack, x_exact)
            
            # Python GTH
            t0 = time.perf_counter()
            Lp, Up = gth_python(A, v.copy(), w.copy())
            t_gth_p += (time.perf_counter() - t0) 
            x_gth_p = solve_lu(Lp, Up, b)
            err_gth_p = rel_error(x_gth_p, x_exact)

        avgt_gth_f = t_gth_f / num_repeats
        avgt_dgetf2 = t_dgetf2 / num_repeats
        avgt_dmmtrf = t_dmmtrf / num_repeats
        avgt_lapack = t_lapack / num_repeats
        avgt_gth_p = t_gth_p / num_repeats
        print(f"{n:6d} "
              f"{avgt_gth_f:12.6f} {err_gth_f:12.2e} "
              f"{avgt_dgetf2:12.6f} {err_dgetf2:12.2e} "
              f"{avgt_dmmtrf:12.6f} {err_dmmtrf:12.2e} "
              f"{avgt_lapack:12.6f} {err_lapack:12.2e} "
              f"{avgt_gth_p:12.6f} {err_gth_p:12.2e}")

if __name__ == "__main__":
    compile_fortran_modules()
    run_test([128,512,1024,2048], num_repeats=5)




 


       



       
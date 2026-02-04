"""
This script tests LU/GTH-style factorizations for NON-SQUARE matrices.

What changes compared to the square case:
-----------------------------------------
1. A is rectangular (m x n).
2. We do NOT solve Ax=b (not defined for rectangular matrices).
3. Accuracy is measured using the LU residual:
       ||A - L U|| / ||A||
4. GTH Python implementation is skipped (square-only algorithm).

Tested routines:
----------------
1. GTH algorithm implemented in Fortran (DMMTF2).
2. Blocked GTH algorithm implemented in Fortran (DMMTRF).
3. LAPACK DGETF2.
4. LAPACK DGETRF.
"""

import numpy as np
import time
import scipy.linalg as la
from scipy.linalg.lapack import get_lapack_funcs
import mpmath as mp
import os
import sys
import subprocess

# ------------------------------------------------------------
# Compile Fortran modules (once)
# ------------------------------------------------------------
def compile_fortran_modules():
    """
    Compiles Fortran sources into Python-callable shared libraries
    using NumPy f2py. This is done only if the .so files do not exist.
    """

    if not os.path.exists("dmmtf2.so"):
        cmd = [sys.executable, "-m", "numpy.f2py", "-c",
               "dmmtf2.f", "-m", "dmmtf2",
               "--fcompiler=gnu95", "-lopenblas"]
        subprocess.run(cmd, check=True)

    if not os.path.exists("dgetf2.so"):
        cmd = [sys.executable, "-m", "numpy.f2py", "-c",
               "dgetf2.f", "-m", "dgetf2",
               "--fcompiler=gnu95", "-lopenblas"]
        subprocess.run(cmd, check=True)

    if not os.path.exists("dmmtrf.so"):
        cmd = [sys.executable, "-m", "numpy.f2py", "-c",
               "dmmtf2.f", "dmmtrf.f",
               "-m", "dmmtrf",
               "--fcompiler=gnu95", "-lopenblas"]
        subprocess.run(cmd, check=True)

# ------------------------------------------------------------
# Rectangular M-matrix generator
# ------------------------------------------------------------
def m_matrix(m, n, seed):
    """
    Generates a rectangular M-matrix-like object.

    Notes:
    ------
    - Off-diagonal entries are non-positive
    - Diagonal dominance is enforced on the leading min(m,n) block
    - For rectangular matrices, the M-matrix concept is weaker,
      but this construction is sufficient for LU testing.
    """
    rng = np.random.default_rng(seed)
    A = -rng.uniform(0, 1e-1, (m, n))

    for i in range(min(m, n)):
        A[i, i] = np.sum(np.abs(A[i, :])) + 1.0

    return A

# ------------------------------------------------------------
# Exact RHS construction: w = A v
# ------------------------------------------------------------
def exact_rhs(A):
    """
    Constructs a nonnegative vector v and computes w = A v exactly
    using high-precision arithmetic.

    Dimensions:
    -----------
    A : (m x n)
    v : (n)
    w : (m)
    """
    m, n = A.shape
    mp.mp.dps = 100

    v_mp = mp.ones(n, 1)
    A_mp = mp.matrix(A.tolist())
    w_mp = A_mp * v_mp

    v = np.ones(n, dtype=np.float64)
    w = np.array(w_mp.tolist(), dtype=np.float64).ravel()
    return v, w

# ------------------------------------------------------------
# LU residual for rectangular matrices
# ------------------------------------------------------------
def lu_residual(A, LU):
    """
    Computes the relative Frobenius residual:

        ||A - L U||_F / ||A||_F

    For rectangular LU:
    -------------------
    - L is m x k with unit diagonal
    - U is k x n
    - k = min(m, n)
    """
    m, n = A.shape
    k = min(m, n)

    L = np.tril(LU[:, :k], -1) + np.eye(m, k)
    U = np.triu(LU[:k, :])

    return la.norm(A - L @ U, ord='fro') / la.norm(A, ord='fro')

# ------------------------------------------------------------
# Fortran GTH (DMMTF2)
# ------------------------------------------------------------
def gth_fortran(A, v, w, module):
    m, n = A.shape
    A_f = np.asfortranarray(A.copy())
    info = np.array(0, dtype=np.int32)

    module.dmmtf2(m=m, n=n, a=A_f, lda=m, v=v, w=w, info=info)
    return A_f

# ------------------------------------------------------------
# Fortran blocked GTH (DMMTRF)
# ------------------------------------------------------------
def dmmtrf_fortran(A, v, w, module, nb):
    m, n = A.shape
    A_f = np.asfortranarray(A.copy())
    info = np.array(0, dtype=np.int32)
    nb_val = np.array(nb, dtype=np.int32)
    work = np.zeros(nb, dtype=np.float64, order="F")

    module.dmmtrf(m=m, n=n, a=A_f, lda=m,
                  v=v, w=w, work=work, nb=nb_val, info=info)
    return A_f

# ------------------------------------------------------------
# LAPACK DGETF2 wrapper
# ------------------------------------------------------------
def dgetf2_fortran(A, module):
    m, n = A.shape
    A_f = np.asfortranarray(A.copy())
    ipiv = np.zeros(min(m, n), dtype=np.int32)
    info = np.array(0, dtype=np.int32)

    module.dgetf2(m=m, n=n, a=A_f, lda=m, ipiv=ipiv, info=info)
    return A_f

# ------------------------------------------------------------
# Main rectangular test driver
# ------------------------------------------------------------
def run_test(sizes, num_repeats):
    import dmmtf2
    import dmmtrf
    import dgetf2

    print(f"{'m':>6} {'n':>6} "
          f"{'dmmtf2 res':>14} "
          f"{'dmmtrf res':>14} "
          f"{'dgetf2 res':>14} "
          f"{'dgetrf res':>14}")
    print("-" * 80)

    for m, n in sizes:
        A = m_matrix(m, n, seed=1)
        v, w = exact_rhs(A)
        nb = 128

        res_dmmtf2 = 0.0
        res_dmmtrf = 0.0
        res_dgetf2 = 0.0
        res_dgetrf = 0.0

        for _ in range(num_repeats):
            # --- GTH unblocked ---
            LU = gth_fortran(A, v.copy(), w.copy(), dmmtf2)
            res_dmmtf2 += lu_residual(A, LU)

            # --- GTH blocked ---
            LU = dmmtrf_fortran(A, v.copy(), w.copy(), dmmtrf, nb)
            res_dmmtrf += lu_residual(A, LU)

            # --- LAPACK DGETF2 ---
            LU = dgetf2_fortran(A, dgetf2)
            res_dgetf2 += lu_residual(A, LU)

            # --- LAPACK DGETRF ---
            A_f = np.asfortranarray(A.copy())
            dgetrf = get_lapack_funcs("getrf", (A_f,))
            LU, piv, _ = dgetrf(A_f)
            res_dgetrf += lu_residual(A, LU)

        print(f"{m:6d} {n:6d} "
              f"{res_dmmtf2/num_repeats:14.2e} "
              f"{res_dmmtrf/num_repeats:14.2e} "
              f"{res_dgetf2/num_repeats:14.2e} "
              f"{res_dgetrf/num_repeats:14.2e}")

# ------------------------------------------------------------
# Entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    compile_fortran_modules()

    # Rectangular test cases (tall matrices)
    test_sizes = [  (512, 256), (1024, 512)  ]

    run_test(test_sizes, num_repeats=1)

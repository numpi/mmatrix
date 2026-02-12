"""
Python GTH wrapper
"""

import subprocess
import sys

import numpy as np

def compile_fortran_modules():
    cmd = [sys.executable, "-m", "numpy.f2py", "-c",
           "-m", "mmatrixlapack", "--fcompiler=gnu95", "-lopenblas",
           "dmmtf2.f", "dmmtrf2.f", "dmmtrf.f", "dmmtrff.f", "dgetf2.f", "dgetrf2.f"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:    
        print("Fortran code compilation failed:")
        print(res.stdout)
        print(res.stderr)
        raise RuntimeError

try:
    import mmatrixlapack
except ImportError:
    compile_fortran_modules()
    import mmatrixlapack

# ------------------------------------------------------------
# GTH Fortran  
def gth_fortran(A, v, w):
    m, n = A.shape
    A_f = np.asfortranarray(A)
    info = np.array(0, dtype=np.int32)
    mmatrixlapack.dmmtf2(m=m, n=n, a=A_f, lda=m, v=v, w=w, info=info)
    return A_f, info[()]

# ------------------------------------------------------------
# DGETF2 Fortran  
def dgetf2_fortran(A):
    m, n = A.shape
    A_f = np.asfortranarray(A.copy())
    ipiv = np.zeros(min(m, n), dtype=np.int32)
    info = np.array(0, dtype=np.int32)
    mmatrixlapack.dgetf2(m=m, n=n, a=A_f, lda=m, ipiv=ipiv, info=info)
    return A_f, ipiv, info[()]

# ------------------------------------------------------------
# DMMTRF Fortran via dmmtf2
def dmmtrf_fortran(A, v, w, nb):
    m, n = A.shape
    A_f = np.asfortranarray(A)
    info = np.array(0, dtype=np.int32)
    nb_val = np.array(nb, dtype=np.int32)
    work = np.zeros(nb, dtype=np.float64, order="F")
    mmatrixlapack.dmmtrf(m=m, n=n, a=A_f, lda=m, v=v, w=w, work=work, nb=nb_val, info=info)
    return A_f, info[()]

# ------------------------------------------------------------
# DMMTRF Fortran via dmmtrf2
def dmmtrff_fortran(A, v, w, nb):
    m, n = A.shape
    A_f = np.asfortranarray(A.copy())
    info = np.array(0, dtype=np.int32)
    nb_val = np.array(nb, dtype=np.int32)
    work_size = min(m, n)
    work = np.zeros(work_size, dtype=np.float64, order="F")
    mmatrixlapack.dmmtrff(m=m, n=n, a=A_f, lda=m, v=v, w=w, work=work, nb=nb_val, info=info)
    return A_f, info[()]

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

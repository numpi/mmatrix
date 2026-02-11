"""
This script tests and compares LU and GTH algorithms 
for FEM matrices with varying condition numbers.

ALGORITHMS COMPARED:
====================
DMMTF2: Unblocked version of the GTH algorithm 
DMMTRF: Blocked version of the GTH algorithm with DMMTF2 inside
DMMTRFF: Blocked version of the GTH algorithm with DMMTRF2 inside
DGETF2: LAPACK unblocked LU factorization
DGETRF: LAPACK blocked LU factorization with DGETRF2 inside
GTH-P: Pure Python implementation of GTH
"""
import numpy as np
import time
import mpmath as mp
import scipy.linalg as la
from scipy.linalg.lapack import get_lapack_funcs
import sys
import subprocess
nb = 16  # block size
n = 64  # matrix size
num_sigmas = 12 # number of different sigma values to test (varying condition numbers)
mp.mp.dps = 32 # set decimal places for mpmath to 32
# ------------------------------------------------------------
# Generate FEM matrix 
def generate_fem_matrix(n, sigma, gamma=0):
    h = 1.0 / (n + 1)
    # Create tridiagonal matrix
    A = np.diag(2 * np.ones(n)) + np.diag(-1 * np.ones(n-1), 1) + np.diag(-1 * np.ones(n-1), -1)
    # Modify first and last diagonal entries
    A[0, 0] = 1 + h * sigma
    A[-1, -1] = 1 + h * sigma
    # Create w 
    w = np.zeros(n)
    w[0] = h * sigma
    w[-1] = h * sigma
    # Create right-hand side b
    grid = h * np.arange(1, n+1)
    b = h * np.sin(np.pi * grid)
    b[0] = b[0] / 2 - h * gamma
    b[-1] = b[-1] / 2 + h * gamma
    # v is ones vector for GTH
    v = np.ones(n)
    return A, v, w, b
# ------------------------------------------------------------
# Compile Fortran modules (we can change setting to compile once)
def compile_fortran_modules():
    cmd = [sys.executable, "-m", "numpy.f2py", "-c", "dmmtrff.f",
           "-m", "dmmtrff", "--fcompiler=gnu95", "-lopenblas"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print("dmmtrff compilation failed:")
        print(res.stderr)   
    cmd = [sys.executable, "-m", "numpy.f2py", "-c", "dmmtf2.f",
           "-m", "dmmtf2", "--fcompiler=gnu95", "-lopenblas"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print("dmmtf2 compilation failed:")
        print(res.stderr)
    cmd = [sys.executable, "-m", "numpy.f2py", "-c", "dgetf2.f",
           "-m", "dgetf2", "--fcompiler=gnu95", "-lopenblas"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print("dgetf2 compilation failed:")
        print(res.stderr)
    cmd = [sys.executable, "-m", "numpy.f2py", "-c", 
           "dmmtf2.f", "dmmtrf.f",  # Include both source files
           "-m", "dmmtrf", 
           "--fcompiler=gnu95", 
           "-lopenblas"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print("dmmtrf compilation failed:")
        print(res.stderr)
              
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
# DMMTRF Fortran via dmmtf2
def dmmtrf_fortran(A, v, w, module, nb):
    m, n = A.shape
    A_f = np.asfortranarray(A)
    info = np.array(0, dtype=np.int32)
    nb_val = np.array(nb, dtype=np.int32)
    work = np.zeros(nb, dtype=np.float64, order="F")
    module.dmmtrf(m=m, n=n, a=A_f, lda=m, v=v, w=w, work=work, nb=nb_val, info=info)
    return A_f

# ------------------------------------------------------------
# DMMTRF Fortran via dmmtrf2
def dmmtrff_fortran(A, v, w, module, nb):
    m, n = A.shape
    A_f = np.asfortranarray(A.copy())
    info = np.array(0, dtype=np.int32)
    nb_val = np.array(nb, dtype=np.int32)
    work_size = min(m, n)
    work = np.zeros(work_size, dtype=np.float64, order="F")
    module.dmmtrff(m=m, n=n, a=A_f, lda=m, v=v, w=w, work=work, nb=nb_val, info=info)
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

# ------------------------------------------------------------
# Solve LUx = b
def solve_lu(L, U, b):
    y = la.solve_triangular(L, b, lower=True, check_finite=False)
    return la.solve_triangular(U, y, check_finite=False)

# ------------------------------------------------------------
# Relative error
def rel_error(x, x_exact):
    return np.max(np.abs(x - x_exact) / np.abs(x_exact))

# ------------------------------------------------------------
# Compute Skeel condition number
def skeel_condition_number(A, x_exact):
    n = A.shape[0]
    # Compute |inv(A)| * |A| * |x|
    abs_A = np.abs(A)
    abs_x = np.abs(x_exact).reshape(-1, 1)
    inv_A = np.linalg.inv(A)
    abs_inv_A = np.abs(inv_A)
    # Compute numerator: || |inv(A)| * |A| * |x| ||_∞
    numerator = np.linalg.norm(abs_inv_A @ (abs_A @ abs_x), np.inf)
    # Compute denominator: ||x||_∞
    denominator = np.linalg.norm(x_exact, np.inf)
    return float(numerator / denominator)

# ------------------------------------------------------------
# Main test with varying condition numbers
def run_test_varying_condition(n, num_sigmas):
    import dmmtf2
    import dgetf2
    import dmmtrf
    import dmmtrff
    
    # Test different sigma values (varying condition numbers)
    sigmas = 10.0 ** (-np.arange(1, num_sigmas + 1))
    
    # Initialize arrays for results
    cond_numbers = []
    skeel_conds = []
    errors_backslash = []
    errors_dmmtf2 = []
    errors_dgetf2 = []
    errors_dmmtrf = []
    errors_dmmtrff = []
    errors_dgetrf = []
    errors_gth_p = []

    
    print("\nTesting FEM matrices with varying condition numbers")
    print("=" * 140)
    print(f"{'Sigma':>10} {'cond(A)':>12} {'skeel(A,x)':>12} "
          f"{'backslash':>12} {'time(s)':>8} "
          f"{'dmmtf2':>12} {'time(s)':>8} "
          f"{'dgetf2':>12} {'time(s)':>8} "
          f"{'dmmtrf':>12} {'time(s)':>8} "
          f"{'dmmtrff':>12} {'time(s)':>8} "
          f"{'dgetrf':>12} {'time(s)':>8} "
          f"{'GTH-P':>12} {'time(s)':>8}")
    print("-" * 140)
    
    for k, sigma in enumerate(sigmas):
        # Generate matrix and vectors
        A, v, w, b = generate_fem_matrix(n, sigma)
        
       
        A_mp = mp.matrix(A.tolist())
        v_mp = mp.matrix(v.reshape(-1, 1).tolist())  
        w_mp = mp.matrix(w.reshape(-1, 1).tolist())   
        b_mp = mp.matrix(b.reshape(-1, 1).tolist())
        
        Avpa = A_mp.copy()
        for i in range(n):
            Avpa[i, i] = 0
        
        f = w_mp - Avpa * v_mp
        for i in range(n):
            Avpa[i, i] = f[i] / v_mp[i]
        x_vpa = np.array((Avpa ** -1 * b_mp).tolist()).ravel()
        
        
        # Compute condition numbers 
        cond_A = float(np.linalg.cond(A))
        skeel_cond = float(skeel_condition_number(A, x_vpa))
        
        cond_numbers.append(cond_A)
        skeel_conds.append(skeel_cond)
        
        # 1. numpy.linalg.solve (backslash equivalent)
        t0 = time.perf_counter()
        x_backslash = np.linalg.solve(A, b)
        t_backslash = time.perf_counter() - t0
        err_backslash = float(rel_error(x_backslash, x_vpa))
        
        # 2. DMMTF2 (unblocked GTH)
        t0 = time.perf_counter()
        A_f = gth_fortran(A, v.copy(), w.copy(), dmmtf2)
        t_dmmtf2 = time.perf_counter() - t0
        Lf = np.tril(A_f, -1) + np.eye(n)
        Uf = np.triu(A_f)
        x_dmmtf2 = solve_lu(Lf, Uf, b)
        err_dmmtf2 = float(rel_error(x_dmmtf2, x_vpa))
        
        # 3. DGETF2 (unblocked LU)
        t0 = time.perf_counter()
        lu_f2, ipiv_f2, _ = dgetf2_fortran(A, dgetf2)
        t_dgetf2 = time.perf_counter() - t0
        L = np.tril(lu_f2, -1) + np.eye(n)
        U = np.triu(lu_f2)
        b_perm = b.copy()
        for i in range(len(ipiv_f2)):
            if ipiv_f2[i] != i + 1:  
                b_perm[[i, ipiv_f2[i]-1]] = b_perm[[ipiv_f2[i]-1, i]]
        y = la.solve_triangular(L, b_perm, lower=True, check_finite=False)
        x_dgetf2 = la.solve_triangular(U, y, check_finite=False)
        err_dgetf2 = float(rel_error(x_dgetf2, x_vpa))
        

        # 4. DMMTRF (blocked GTH using dmmtf2)
        t0 = time.perf_counter()
        A_trfa = dmmtrf_fortran(A, v.copy(), w.copy(), dmmtrf, nb)
        t_dmmtrf = time.perf_counter() - t0
        L = np.tril(A_trfa, -1) + np.eye(n)
        U = np.triu(A_trfa)
        x_dmmtrf = solve_lu(L, U, b)
        err_dmmtrf = float(rel_error(x_dmmtrf, x_vpa))
        
        # 5. DMMTRFF (blocked GTH using dmmtrf2)
        t0 = time.perf_counter()
        A_trff, info = dmmtrff_fortran(A, v.copy(), w.copy(), dmmtrff, nb)
        t_dmmtrff = time.perf_counter() - t0
        L = np.tril(A_trff, -1) + np.eye(n)
        U = np.triu(A_trff)
        x_dmmtrff = solve_lu(L, U, b)
        err_dmmtrff = float(rel_error(x_dmmtrff, x_vpa))
        
        # 6. DGETRF (blocked LU)
        A_f = np.asfortranarray(A)
        dgetrf = get_lapack_funcs('getrf', (A_f,))
        
        t0 = time.perf_counter()
        lu_rf, piv_rf, _ = dgetrf(A_f)
        t_dgetrf = time.perf_counter() - t0
        x_dgetrf = la.lu_solve((lu_rf, piv_rf), b)
        err_dgetrf = float(rel_error(x_dgetrf, x_vpa))
        
        # 7. Python GTH
        t0 = time.perf_counter()
        Lp, Up = gth_python(A, v.copy(), w.copy())
        t_gth_p = time.perf_counter() - t0
        x_gth_p = solve_lu(Lp, Up, b)
        err_gth_p = float(rel_error(x_gth_p, x_vpa))
        
        # Store errors and times
        errors_backslash.append(err_backslash)
        errors_dmmtf2.append(err_dmmtf2)
        errors_dgetf2.append(err_dgetf2)
        errors_dmmtrf.append(err_dmmtrf)
        errors_dmmtrff.append(err_dmmtrff)
        errors_dgetrf.append(err_dgetrf)
        errors_gth_p.append(err_gth_p)
    
     
        print(f"{sigma:10.1e} {cond_A:12.2e} {skeel_cond:12.2e} "
              f"{err_backslash:12.2e} {t_backslash*1000:8.4f} "
              f"{err_dmmtf2:12.2e} {t_dmmtf2*1000:8.4f} "
              f"{err_dgetf2:12.2e} {t_dgetf2*1000:8.4f} "
              f"{err_dmmtrf:12.2e} {t_dmmtrf*1000:8.4f} "
              f"{err_dmmtrff:12.2e} {t_dmmtrff*1000:8.4f} "
              f"{err_dgetrf:12.2e} {t_dgetrf*1000:8.4f} "
              f"{err_gth_p:12.2e} {t_gth_p*1000:8.4f}")
    
    create_plot(sigmas, cond_numbers, skeel_conds, errors_backslash, errors_dmmtf2, 
                errors_dgetf2, errors_dmmtrf, errors_dmmtrff, errors_dgetrf, 
                errors_gth_p, n)

    
    return {
        'sigmas': sigmas,
        'cond_numbers': cond_numbers,
        'skeel_conds': skeel_conds,
        'errors': {
            'backslash': errors_backslash,
            'dmmtf2': errors_dmmtf2,
            'dgetf2': errors_dgetf2,
            'dmmtrf': errors_dmmtrf,
            'dmmtrff': errors_dmmtrff,
            'dgetrf': errors_dgetrf,
            'gth_p': errors_gth_p
        }
    }

# ------------------------------------------------------------
def create_plot(sigmas, cond_numbers, skeel_conds, errors_backslash, errors_dmmtf2,
                errors_dgetf2, errors_dmmtrf, errors_dmmtrff, errors_dgetrf,
                errors_gth_p, n):
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(12, 8))
    
    indices = np.arange(1, len(sigmas) + 1)
    
    plt.semilogy(indices, np.array(cond_numbers) * np.finfo(float).eps, '--', 
                    linewidth=2, markersize=8, label='κ(A)·u')
    
    plt.semilogy(indices, np.array(skeel_conds) * np.finfo(float).eps, 'xr', 
                    linewidth=1, markersize=6, label='skeel cond(A,x)')
    
    
    plt.semilogy(indices, errors_backslash, 'sb', 
                    linewidth=1, markersize=10, label='relerr-norm-backslash')
    
    plt.semilogy(indices, errors_dmmtf2, 'dc', 
                    linewidth=2, markersize=8, label='relerr-norm-dmmtf2 (GTH)')
    
    plt.semilogy(indices, errors_dgetf2, '^g', 
                    linewidth=1, markersize=8, label='relerr-norm-dgetf2')
    
    plt.semilogy(indices, errors_dmmtrf, 'vm', 
                    linewidth=2, markersize=8, label='relerr-norm-dmmtrf (blocked GTH via dmmtrf2)')
    
    plt.semilogy(indices, errors_dmmtrff, '*y', 
                    linewidth=2, markersize=10, label='relerr-norm-dmmtrff (blocked GTH via dmmtrf2)')
    
    plt.semilogy(indices, errors_dgetrf, '.r', 
                    linewidth=4, markersize=8, label='relerr-norm-dgetrf')
    
    plt.semilogy(indices, errors_gth_p, 'ok', 
                    linewidth=1, markersize=6, label='relerr-norm-GTH-Python')
    
    plt.xlabel('Index (decreasing σ)', fontsize=12)
    plt.ylabel('Relative Error / Condition Bound', fontsize=12)
    plt.title(f'{n}x{n} linear systems with FEM matrices, varying condition number', 
                fontsize=14, fontweight='bold')
    
    plt.legend(loc='best', fontsize=10)
    
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    plt.savefig(f'fem_matrix_comparison_n{n}.png', dpi=300, bbox_inches='tight')
    print(f"\nPlot saved as 'fem_matrix_comparison_n{n}.png'")
    
    plt.show()

# ------------------------------------------------------------
if __name__ == "__main__":
    compile_fortran_modules()
    results = run_test_varying_condition(n, num_sigmas)
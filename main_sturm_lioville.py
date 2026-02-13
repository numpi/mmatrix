"""
1D FEM to solve SL problem with Robin BC in FEniCSx
Compares LU and GTH algorithms using Fortran implementations
"""
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import pandas as pd
import time
import sys
import subprocess
import os
from dolfinx import mesh, fem
from dolfinx.fem import petsc
from dolfinx.fem.petsc import LinearProblem
from ufl import TrialFunction, TestFunction, grad, inner, dx, ds
from scipy.linalg import lu_factor, lu_solve

from mmatrix import *

# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------
a_val, b_val, sigma = 0, 1, 1e-9
k_a = k_b = sigma
g_a = g_b = 0

# Test with increasing N
Ns = [256,512,1024,2048,4096]
nb=128

# ------------------------------------------------------------
# Symbolic exact solution
# ------------------------------------------------------------
x_sym = sp.symbols('x')
p_sym = 1 + x_sym
q_sym = x_sym / 100

A_sym, B_sym = sp.symbols('A B')
u_sym = (1 - x_sym)**2 * sp.sin(sp.pi*x_sym)/(sp.pi)**2 + A_sym*x_sym + B_sym
up_sym = sp.diff(u_sym, x_sym)

# Solve Robin BCs
u0 = u_sym.subs(x_sym, 0)
up0 = up_sym.subs(x_sym, 0)
u1 = u_sym.subs(x_sym, 1)
up1 = up_sym.subs(x_sym, 1)

eq1 = sp.Eq(-p_sym.subs(x_sym,0)*up0 + k_a*u0, g_a)
eq2 = sp.Eq(p_sym.subs(x_sym,1)*up1 + k_b*u1, g_b)

solAB = sp.solve([eq1, eq2], [A_sym, B_sym], dict=True)[0]
A_val = float(solAB[A_sym])
B_val = float(solAB[B_sym])

# Compute f(x) and exact solution
f_sym = -sp.diff(p_sym*up_sym, x_sym) + q_sym*u_sym
f_sub = f_sym.subs({A_sym: A_val, B_sym: B_val})
u_ex_sym = u_sym.subs({A_sym: A_val, B_sym: B_val})

# Convert to numpy functions
f_np = sp.lambdify(x_sym, f_sub, 'numpy')
u_np = sp.lambdify(x_sym, u_ex_sym, 'numpy')
p_np = sp.lambdify(x_sym, p_sym, 'numpy')
q_np = sp.lambdify(x_sym, q_sym, 'numpy')

# ------------------------------------------------------------
def solve_lu(L, U, b, ipiv=None):
    """Solve LUx = b"""
    n = L.shape[0]
    if ipiv is not None:
        b_perm = b.copy()
        for i in range(len(ipiv)):
            if ipiv[i] != i + 1:
                j = ipiv[i] - 1
                b_perm[i], b_perm[j] = b_perm[j], b_perm[i]
        y = np.zeros(n, dtype=np.float64)
        for i in range(n):
            y[i] = b_perm[i] - np.dot(L[i, :i], y[:i])
    else:
        y = np.zeros(n, dtype=np.float64)
        for i in range(n):
            y[i] = b[i] - np.dot(L[i, :i], y[:i])
    
    x = np.zeros(n, dtype=np.float64)
    for i in range(n-1, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i+1:], x[i+1:])) / U[i, i]
    return x
 
# Store results
results = {
    'N': [], 'h': [], 'cond': [],
    'LU_FEM_time': [], 'LU_FEM_error': [],
    'DMMTF2_time': [], 'DMMTF2_error': [],
    'DMMTRF_time': [], 'DMMTRF_error': [],
    'DMMTRFF_time': [], 'DMMTRFF_error': [],
    'DGETF2_time': [], 'DGETF2_error': [],
    'DGETRF_time': [], 'DGETRF_error': []
}

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------
# Compile Fortran modules
if not compile_fortran_modules():
    print("Failed to compile Fortran modules. Exiting.")
    sys.exit(1)

print("\nFortran modules compiled successfully. Starting computations...\n")

for N in Ns: 
    try:
        # Create mesh
        domain = mesh.create_interval(MPI.COMM_WORLD, N, [a_val, b_val])
        V = fem.functionspace(domain, ("Lagrange", 1))
        
        # Boundary setup
        tdim = domain.topology.dim
        fdim = tdim - 1
        domain.topology.create_connectivity(fdim, tdim)
        
        boundary_facets = mesh.locate_entities_boundary(domain, fdim,
                           lambda x: np.isclose(x[0], a_val) | np.isclose(x[0], b_val))
        
        marked_facets, marked_values = [], []
        for facet in boundary_facets:
            x = domain.geometry.x[facet]
            if np.isclose(x[0], a_val):
                marked_facets.append(facet)
                marked_values.append(1)
            else:
                marked_facets.append(facet)
                marked_values.append(2)
        
        marked_facets = np.array(marked_facets, dtype=np.int32)
        marked_values = np.array(marked_values, dtype=np.int32)
        boundary_tags = mesh.meshtags(domain, fdim, marked_facets, marked_values)
        ds = ds(domain, subdomain_data=boundary_tags)
        
        # Functions
        V2 = fem.functionspace(domain, ("Lagrange", 2))
        uex = fem.Function(V2, name="u_exact")
        u_ex = fem.Function(V)
        fx = fem.Function(V)
        px = fem.Function(V)
        qx = fem.Function(V)
        
        uex.interpolate(lambda x: u_np(x[0]))
        u_ex.interpolate(lambda x: u_np(x[0]))
        fx.interpolate(lambda x: f_np(x[0]))
        px.interpolate(lambda x: p_np(x[0]))
        qx.interpolate(lambda x: q_np(x[0]))
        
        # Variational forms
        u = TrialFunction(V)
        v = TestFunction(V)
        
        a = (inner(px*grad(u), grad(v))*dx + inner(qx*u, v)*dx + 
             k_a*u*v*ds(1) + k_b*u*v*ds(2))
        L = (inner(fx, v)*dx + g_a*v*ds(1) + g_b*v*ds(2))
        
        # Assemble system
        A = petsc.create_matrix(fem.form(a))
        petsc.assemble_matrix(A, fem.form(a))
        A.assemble()
        
        b = petsc.create_vector(fem.form(L))
        petsc.assemble_vector(b, fem.form(L))
        b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        b.assemble()
        
        # Assemble w_gth
        w_form = fem.form(qx*u*dx + k_a*u*ds(1) + k_b*u*ds(2))
        w_gth = petsc.create_vector(w_form)
        petsc.assemble_vector(w_gth, w_form)
        w_gth.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        w_gth.assemble()
        
        # Convert to numpy
        A_np = A.convert("dense").getDenseArray().astype(np.float64)
        b_np = b.array.astype(np.float64).copy()
        w_gth_np = w_gth.array.astype(np.float64).copy()
        
        # Condition number
        cond_num = np.linalg.cond(A_np)
        h_val = (b_val - a_val) / N
        v_vec = np.ones(A_np.shape[0], dtype=np.float64)
        
        # Store basic info
        results['N'].append(N)
        results['h'].append(h_val)
        results['cond'].append(cond_num)
        
        # L2 norm of exact solution
        L2_exact = fem.form(inner(uex, uex)*dx)
        error_locale = fem.assemble_scalar(L2_exact)
        exact_L2 = np.sqrt(domain.comm.allreduce(error_locale, op=MPI.SUM))
        
        # ----------------------------------------------------
        # Method 1: LU via FEniCSx
        # ----------------------------------------------------
        t0 = time.perf_counter()
        problem = LinearProblem(a, L, petsc_options={
            "ksp_type": "preonly", 
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps"
        })
        uh_lu = problem.solve()
        t_lu_fem = time.perf_counter() - t0
        results['LU_FEM_time'].append(t_lu_fem)
        
        # Compute error
        L2_lu_fem = fem.form(inner(uh_lu-uex, uh_lu-uex)*dx)
        error_local = fem.assemble_scalar(L2_lu_fem)
        error_lu_fem = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM)) / exact_L2
        results['LU_FEM_error'].append(error_lu_fem)
        
        
        # ----------------------------------------------------
        # Method 3: DMMTF2 (Fortran)
        # ----------------------------------------------------
        t0 = time.perf_counter()
        A_dmmtf2, _ = dmmtf2_fortran(A_np, v_vec, w_gth_np.copy())
        t_dmmtf2 = time.perf_counter() - t0
        L_dmmtf2 = np.tril(A_dmmtf2, -1) + np.eye(A_np.shape[0])
        U_dmmtf2 = np.triu(A_dmmtf2)
        x_dmmtf2 = solve_lu(L_dmmtf2, U_dmmtf2, b_np)

        results['DMMTF2_time'].append(t_dmmtf2)
        
        uh_dmmtf2 = fem.Function(V)
        uh_dmmtf2.x.array[:] = x_dmmtf2
        
        L2_dmmtf2 = fem.form(inner(uh_dmmtf2-uex, uh_dmmtf2-uex)*dx)
        error_local = fem.assemble_scalar(L2_dmmtf2)
        error_dmmtf2 = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM)) / exact_L2
        results['DMMTF2_error'].append(error_dmmtf2)
        
        # ----------------------------------------------------
        # Method 4: DMMTRF (Fortran)
        # ----------------------------------------------------
        t0 = time.perf_counter()
        A_dmmtrf,_ = dmmtrf_fortran(A_np, v_vec, w_gth_np.copy(), nb)
        t_dmmtrf = time.perf_counter() - t0
        L_dmmtrf = np.tril(A_dmmtrf, -1) + np.eye(A_np.shape[0])
        U_dmmtrf = np.triu(A_dmmtrf)
        x_dmmtrf = solve_lu(L_dmmtrf, U_dmmtrf, b_np)
        results['DMMTRF_time'].append(t_dmmtrf)
        
        uh_dmmtrf = fem.Function(V)
        uh_dmmtrf.x.array[:] = x_dmmtrf
        
        L2_dmmtrf = fem.form(inner(uh_dmmtrf-uex, uh_dmmtrf-uex)*dx)
        error_local = fem.assemble_scalar(L2_dmmtrf)
        error_dmmtrf = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM)) / exact_L2
        results['DMMTRF_error'].append(error_dmmtrf)
        # ----------------------------------------------------
        # Method 5: DMMTRFF (Fortran)
        # ----------------------------------------------------
        t0 = time.perf_counter()
        A_dmmtrff, _ = dmmtrff_fortran(A_np, v_vec, w_gth_np.copy(), nb)
        t_dmmtrff = time.perf_counter() - t0
        L_dmmtrff = np.tril(A_dmmtrff, -1) + np.eye(A_np.shape[0])
        U_dmmtrff = np.triu(A_dmmtrff)
        x_dmmtrff = solve_lu(L_dmmtrff, U_dmmtrff, b_np)

        results['DMMTRFF_time'].append(t_dmmtrff)
        
        uh_dmmtrff = fem.Function(V)
        uh_dmmtrff.x.array[:] = x_dmmtrff
        
        L2_dmmtrff = fem.form(inner(uh_dmmtrff-uex, uh_dmmtrff-uex)*dx)
        error_local = fem.assemble_scalar(L2_dmmtrff)
        error_dmmtrff = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM)) / exact_L2
        results['DMMTRFF_error'].append(error_dmmtrff)
        
        # ----------------------------------------------------
        # Method 6: DGETF2 (Fortran)
        # ----------------------------------------------------
        t0 = time.perf_counter()
        A_dgetf2, ipiv, _ = dgetf2_fortran(A_np)
        t_dgetf2 = time.perf_counter() - t0
        L_dgetf2 = np.tril(A_dgetf2, -1) + np.eye(A_np.shape[0])
        U_dgetf2 = np.triu(A_dgetf2)
        x_dgetf2 = solve_lu(L_dgetf2, U_dgetf2, b_np, ipiv)
        results['DGETF2_time'].append(t_dgetf2)
            
        uh_dgetf2 = fem.Function(V)
        uh_dgetf2.x.array[:] = x_dgetf2
            
        L2_dgetf2 = fem.form(inner(uh_dgetf2-uex, uh_dgetf2-uex)*dx)
        error_local = fem.assemble_scalar(L2_dgetf2)
        error_dgetf2 = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM)) / exact_L2
        results['DGETF2_error'].append(error_dgetf2)
        
        # ----------------------------------------------------
        # Method 7: DGETRF (LAPACK)
        # ----------------------------------------------------
        t0 = time.perf_counter()
        lu_rf, piv_rf = lu_factor(A_np, check_finite=False)
        t_dgetrf = time.perf_counter() - t0
        x_dgetrf = lu_solve((lu_rf, piv_rf), b_np, check_finite=False)
        results['DGETRF_time'].append(t_dgetrf)
        
        uh_dgetrf = fem.Function(V)
        uh_dgetrf.x.array[:] = x_dgetrf
        
        L2_dgetrf = fem.form(inner(uh_dgetrf-uex, uh_dgetrf-uex)*dx)
        error_local = fem.assemble_scalar(L2_dgetrf)
        error_dgetrf = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM)) / exact_L2
        results['DGETRF_error'].append(error_dgetrf)
        
        # ----------------------------------------------------
        # Print current results
        # ----------------------------------------------------
        print(f"RESULTS for N={N}:")
        print(f"{'-'*50}")
        #print(f"  Condition number: {cond_num:.2e}")
        print(f"  Mesh size h: {h_val:.2e}")
        print(f"\n  {'Method':<15} {'L2 Error':<15} {'Time (ms)':<12}")
        print(f"  {'-'*54}")
        
        print(f"  {'LU-FEM':<15} {error_lu_fem:<15.2e} {t_lu_fem*1000:8.2f}")
        
        if not np.isnan(results['DMMTF2_error'][-1]):
            print(f"  {'DMMTF2':<15} {error_dmmtf2:<15.2e} {t_dmmtf2*1000:8.2f}")
        if not np.isnan(results['DMMTRF_error'][-1]):
            print(f"  {'DMMTRF':<15} {error_dmmtrf:<15.2e} {t_dmmtrf*1000:8.2f}")
        if not np.isnan(results['DMMTRFF_error'][-1]):
            print(f"  {'DMMTRFF':<15} {error_dmmtrff:<15.2e} {t_dmmtrff*1000:8.2f}")
        if not np.isnan(results['DGETF2_error'][-1]):
            print(f"  {'DGETF2':<15} {error_dgetf2:<15.2e} {t_dgetf2*1000:8.2f}")
        if not np.isnan(results['DGETRF_error'][-1]):
            print(f"  {'DGETRF':<15} {error_dgetrf:<15.2e} {t_dgetrf*1000:8.2f}")
        
    except Exception as e:
        print(f"Error processing N={N}: {e}")
        import traceback
        traceback.print_exc()
        continue

# ------------------------------------------------------------
# Create results table
# ------------------------------------------------------------
print("\n" + "-"*100)

df = pd.DataFrame(results)
print("\n", df.to_string(index=False))

# ------------------------------------------------------------
# Create plots
# ------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
x_coords = V.tabulate_dof_coordinates()[:,0]
sort_idx = np.argsort(x_coords)
ax1.loglog(df['h'], df['LU_FEM_error'], 'ko-', linewidth=2, markersize=8, label='LU-FEM')
ax1.loglog(df['h'], df['DMMTF2_error'], 'b^-', linewidth=2, markersize=8, label='DMMTF2')
ax1.loglog(df['h'], df['DMMTRF_error'], 'gd-', linewidth=2, markersize=8, label='DMMTRF')
ax1.loglog(df['h'], df['DMMTRFF_error'], 'm*-', linewidth=2, markersize=8, label='DMMTRFF')
ax1.loglog(df['h'], df['DGETF2_error'], 'cv-', linewidth=2, markersize=8, label='DGETF2')
ax1.loglog(df['h'], df['DGETRF_error'], 'yp-', linewidth=2, markersize=8, label='DGETRF')
ax1.set_xlabel('h (mesh size)')
ax1.set_ylabel('Relative L2 Error')
ax1.invert_xaxis()  
ax1.set_title('Convergence: L2 Error vs mesh size')
ax1.legend(loc='best', fontsize=8)
ax1.grid(True, alpha=0.3)

# Plot 2: CPU Time
ax2.loglog(df['N'], df['LU_FEM_time'], 'ko-', linewidth=2, markersize=8, label='LU-FEM')
ax2.loglog(df['N'], df['DMMTF2_time'], 'b^-', linewidth=2, markersize=8, label='DMMTF2')
ax2.loglog(df['N'], df['DMMTRF_time'], 'gd-', linewidth=2, markersize=8, label='DMMTRF')
ax2.loglog(df['N'], df['DMMTRFF_time'], 'm*-', linewidth=2, markersize=8, label='DMMTRFF')
ax2.loglog(df['N'], df['DGETF2_time'], 'cv-', linewidth=2, markersize=8, label='DGETF2')
ax2.loglog(df['N'], df['DGETRF_time'], 'yp-', linewidth=2, markersize=8, label='DGETRF')
ax2.set_xlabel('N (elements)')
ax2.set_ylabel('CPU Time (s)')
ax2.set_title('CPU time vs problem size')
ax2.legend(loc='best', fontsize=8)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
    
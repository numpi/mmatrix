function solver_demo()
% This tool allows users to compare different factorization methods
% for solving linear systems Ax = b, where A is an M-matrix and 
% and vectors v and w are chosen such that Av = w.
%
clc; close all;

build_mex(); 

fprintf('============================================\n');
fprintf('    MATRIX SOLVER COMPARISON TOOL\n');
fprintf('============================================\n\n');

% Get user preferences
[method_flags, use_custom] = get_user_preferences();

% Get matrix data
if use_custom
    [A, b, v, w] = get_custom_input();
    run_single_matrix(A, b, v, w, method_flags);
else
    run_dimension_study(method_flags);
end

end

%--------------------------------------------------------------------------
function [method_flags, use_custom] = get_user_preferences()
% Get user preferences for methods and problem type

fprintf('Choose problem type:\n');
fprintf('  1. Use predefined FEM problem\n');
fprintf('  2. Input an M-matrix A and vectors v,w,b\n');
choice = input('Enter your choice (1 or 2): ');

use_custom = (choice == 2);

fprintf('\nAvailable solution methods:\n');
fprintf('  1. GTH(pure MATLAB)\n');
fprintf('  2. LU(dgetrf)\n');
fprintf('  3. GTH(DMMTF2)\n');
fprintf('  4. GTH(DMMTRF)\n');
fprintf('  5. GTH(DMMTRFF)\n');
fprintf('\nYou can run multiple methods (e.g., "1,3,5" or "all")\n');

method_input = input('Enter method numbers to run: ', 's');

% Parse method selection
if strcmpi(method_input, 'all')
    method_flags = struct('gth', true, 'lu', true, 'dmmtf2', true, ...
                          'dmmtrf', true, 'dmmtrff', true);
else
    selected = str2num(method_input);
    method_flags = struct('gth', false, 'lu', false, 'dmmtf2', false, ...
                          'dmmtrf', false, 'dmmtrff', false);
    
    if any(selected == 1)
        method_flags.gth = true;
    end
    if any(selected == 2)
        method_flags.lu = true;
    end
    if any(selected == 3)
        method_flags.dmmtf2 = true;
    end
    if any(selected == 4)
        method_flags.dmmtrf = true;
    end
    if any(selected == 5)
        method_flags.dmmtrff = true;
    end
end

fprintf('\n');
end

%--------------------------------------------------------------------------
function [A, b, v, w] = get_custom_input()
% Get custom matrix and vectors from user

fprintf('\n--- Custom Matrix Input ---\n');

% Get matrix A
n = input('Enter matrix dimension n: ');
fprintf('Enter matrix A (row by row, or as MATLAB expression):\n');
fprintf('Example: [2 -1 0; -1 2 -1; 0 -1 2]\n');
A_input = input('A = ');
if ischar(A_input)
    A = eval(A_input);
else
    A = A_input;
end

% Verify matrix size
if size(A,1) ~= n || size(A,2) ~= n
    error('Matrix dimensions do not match specified n');
end

% Get vector b
fprintf('\nEnter right-hand side vector b (size %d x 1):\n', n);
fprintf('Example: [1; 2; 3]\n');
b_input = input('b = ');
if ischar(b_input)
    b = eval(b_input);
else
    b = b_input;
end

% Get vector v (for GTH and DMMT methods)
fprintf('\nEnter vector v (size %d x 1, press Enter for ones):\n', n);
v_input = input('v = (default: ones) ', 's');
if isempty(v_input)
    v = ones(n,1);
else
    v = eval(v_input);
end

% Get vector w (for GTH and DMMT methods)
fprintf('\nEnter vector w (size %d x 1, press Enter for zeros):\n', n);
w_input = input('w = (default: zeros) ', 's');
if isempty(w_input)
    w = zeros(n,1);
else
    w = eval(w_input);
end

fprintf('\nMatrix input complete.\n\n');
end

%--------------------------------------------------------------------------
function run_single_matrix(A, b, v, w, method_flags)
% Run selected methods on a single matrix and display table

fprintf('\n--- Running Single Matrix Analysis ---\n');
fprintf('Matrix size: %d x %d\n', size(A,1), size(A,2));
fprintf('Condition number: %.2e\n\n', cond(A));

% Compute exact solution
x_exact = MatrixSolverFunctions.exact_solution(A, v, w, b);

% Initialize results structure
results = struct();
results.methods = {};
results.errors = [];
results.times = [];

% Run selected methods
if method_flags.gth
    fprintf('Running GTH...\n');
    [time, error] = MatrixSolverFunctions.run_gthlu(A, b, x_exact, v, w);
    results.methods{end+1} = 'GTH';
    results.errors(end+1) = error;
    results.times(end+1) = time;
end

if method_flags.lu
    fprintf('Running Standard LU...\n');
    [time, error] = MatrixSolverFunctions.run_standard_lu(A, b, x_exact);
    results.methods{end+1} = 'Standard LU';
    results.errors(end+1) = error;
    results.times(end+1) = time;
end

if method_flags.dmmtf2
    fprintf('Running DMMTF2...\n');
    [time, error] = MatrixSolverFunctions.run_dmmtf2(A, b, x_exact, v, w);
    results.methods{end+1} = 'DMMTF2';
    results.errors(end+1) = error;
    results.times(end+1) = time;
end

if method_flags.dmmtrf
    nb = input('Enter block size for DMMTRF (default: 64): ');
    if isempty(nb)
        nb = 64;
    end
    fprintf('Running DMMTRF (nb=%d)...\n', nb);
    [time, error] = MatrixSolverFunctions.run_dmmtrf(A, b, x_exact, v, w, nb);
    results.methods{end+1} = 'DMMTRF';
    results.errors(end+1) = error;
    results.times(end+1) = time;
end

if method_flags.dmmtrff
    nb = input('Enter block size for DMMTRFF (default: 64): ');
    if isempty(nb)
        nb = 64;
    end
    fprintf('Running DMMTRFF (nb=%d)...\n', nb);
    [time, error] = MatrixSolverFunctions.run_dmmtrff(A, b, x_exact, v, w, nb);
    results.methods{end+1} = 'DMMTRFF';
    results.errors(end+1) = error;
    results.times(end+1) = time;
end

% Display results table
MatrixSolverFunctions.display_single_results_table(results, cond(A));
end

%--------------------------------------------------------------------------
function run_dimension_study(method_flags)
% Run dimension study with FEM problem and display tables

fprintf('\n--- Running ---\n');

% Get parameters
params = struct();
params.sigma = input('Enter sigma value (default: 1e-3): ');
if isempty(params.sigma)
    params.sigma = 1e-3;
end

ns_input = input('Enter matrix dimensions as array (default: [128 256 512 1024]): ', 's');
if isempty(ns_input)
    params.ns = [128 256 512 1024];
else
    params.ns = eval(ns_input);
end

params.nb = input('Enter block size (default: 64): ');
if isempty(params.nb)
    params.nb = 64;
end

% Run experiments
results = MatrixSolverFunctions.run_dimension_experiments(params, method_flags);

% Display results tables
MatrixSolverFunctions.display_dimension_tables(results, params);
end
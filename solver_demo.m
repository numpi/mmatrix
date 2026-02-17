function solver_demo()
% SOLVER_DEMO compares various factorization methods for solving linear systems
% Ax = b, where A is an M-matrix and vectors v, w satisfy Av = w.
%
% Available methods:
%   • GTH      - Pure MATLAB implementation (uses v,w)
%   • LU       - LAPACK blocked LU factorization (DGETRF) (does not use v,w)
%   • DMMTF2   - Unblocked GTH algorithm (uses v,w)
%   • DMMTRF   - Blocked GTH with DMMTF2 inside (uses v,w)
%   • DMMTRFF  - Blocked GTH with DMMTRF2 inside (uses v,w)
%
% See also: mmatrix, build_mex
%
clc; close all;

build_mex();  % Build MEX files if needed

fprintf('============================================\n');
fprintf('   COMPARE FACTORIZATION METHODS\n');
fprintf('============================================\n\n');

% Get user preferences including problem type and parameters
[method_flags, params, use_custom, use_workspace] = get_user_preferences();

% Get matrix data based on user choice
if use_workspace
    [A, b, v, w] = get_from_workspace();
    run_problem(A, b, v, w, method_flags, params);
elseif use_custom
    % Check if any GTH methods are selected (methods that need v,w)
    gth_methods_selected = method_flags.gth || method_flags.dmmtf2 || ...
                           method_flags.dmmtrf || method_flags.dmmtrff;
    
    [A, b, v, w] = get_custom_input(method_flags, gth_methods_selected);
    run_problem(A, b, v, w, method_flags, params);
else
    run_predefined_test(method_flags, params);
end

end

%--------------------------------------------------------------------------
function [method_flags, params, use_custom, use_workspace] = get_user_preferences()
% Get user preferences for methods and problem type

fprintf('Choose problem type:\n');
fprintf('  1. Use predefined test problem\n');
fprintf('  2. Input an M-matrix A and vectors v, w, b\n');
fprintf('  3. Load A, b, v, w from workspace\n');
choice = input('Enter your choice (1, 2, or 3): ');

use_custom = (choice == 2);
use_workspace = (choice == 3);

[method_flags, params] = get_method_preferences();

fprintf('\n');
end

%--------------------------------------------------------------------------
function [A, b, v, w] = get_from_workspace()
% Load matrix and vectors directly from workspace

fprintf('\n--- Loading from Workspace ---\n');

% Check if variables exist
if ~evalin('base', 'exist(''A'', ''var'')')
    error('Variable ''A'' not found in workspace');
end
if ~evalin('base', 'exist(''b'', ''var'')')
    error('Variable ''b'' not found in workspace');
end

A = evalin('base', 'A');
b = evalin('base', 'b');

% Verify A is a square matrix
if size(A,1) ~= size(A,2)
    error('A must be a square matrix');
end
n = size(A,1);

% Verify vector dimensions
if size(b,1) ~= n
    error('b must have %d rows to match matrix A', n);
end

% Try to load v and w  
if evalin('base', 'exist(''v'', ''var'')')
    v = evalin('base', 'v');
    if size(v,1) ~= n
        error('v must have %d rows to match matrix A', n);
    end
else
    v = ones(n,1);
    fprintf('Warning: v not found in workspace, using ones(%d,1)\n', n);
end

if evalin('base', 'exist(''w'', ''var'')')
    w = evalin('base', 'w');
    if size(w,1) ~= n
        error('w must have %d rows to match matrix A', n);
    end
else
    fprintf('Warning: w not found in workspace\n');
end

end

%--------------------------------------------------------------------------
function [method_flags, params] = get_method_preferences()
% Get user preferences for solution methods and their parameters

fprintf('\nAvailable solution methods:\n');
fprintf('  1. GTH(pure MATLAB) - uses v,w\n');
fprintf('  2. LU(dgetrf) - does not use v,w\n');
fprintf('  3. GTH(DMMTF2) - uses v,w (unblocked)\n');
fprintf('  4. GTH(DMMTRF) - uses v,w (blocked, needs block size)\n');
fprintf('  5. GTH(DMMTRFF) - uses v,w (blocked, needs block size)\n');
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

% Ask for block size if any blocked methods are selected
params = struct();
if method_flags.dmmtrf || method_flags.dmmtrff
    fprintf('\n--- Block Size Configuration ---\n');
    fprintf('Blocked methods (DMMTRF, DMMTRFF) need a block size parameter.\n');
    nb_input = input('Enter block size (default: 64): ', 's');
    if isempty(nb_input)
        params.nb = 64;
    else
        params.nb = str2double(nb_input);
    end
    fprintf('Using block size: %d\n', params.nb);
else
    params.nb = 64;  
end
end

%--------------------------------------------------------------------------
function [A, b, v, w] = get_custom_input(method_flags, gth_methods_selected)
% Get  matrix and vectors from user

fprintf('\n--- Custom Matrix Input ---\n');

% Get matrix A
n = input('Enter matrix dimension n: ');
fprintf('Enter matrix A :\n');
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

% Only ask for v and w if GTH methods are selected
if gth_methods_selected
    % Get  v (for GTH and DMMT methods)
    fprintf('\nEnter vector v (size %d x 1, press Enter for ones):\n', n);
    v_input = input('v = (default: ones) ', 's');
    if isempty(v_input)
        v = ones(n,1);
    else
        v = eval(v_input);
    end

    % Get w (for GTH and DMMT methods)
    fprintf('\nEnter vector w (size %d x 1):\n', n);
    w_input = input('w = ', 's');
    w = eval(w_input);
end

fprintf('\nMatrix input complete.\n\n');
end

%--------------------------------------------------------------------------
function run_problem(A, b, v, w, method_flags, params)
% Run selected methods and display table

fprintf('Matrix size: %d x %d\n', size(A,1), size(A,2));
fprintf('Condition number: %.2e\n\n', cond(A));

% Compute exact solutions
x_exact_lu = mmatrix.exact_solution_lu(A, b);  

% Only compute GTH reference if needed
gth_methods_selected = method_flags.gth || method_flags.dmmtf2 || ...
                       method_flags.dmmtrf || method_flags.dmmtrff;

if gth_methods_selected
    x_exact_gth = mmatrix.exact_solution(A, v, w, b);  
else
    x_exact_gth = [];  
end

 
results = struct();
results.methods = {};
results.errors = [];
results.times = [];
results.solver_type = {};

% Run selected methods
if method_flags.gth
    fprintf('Running GTH (uses v,w)...\n');
    [time, error] = mmatrix.run_gthlu(A, b, x_exact_gth, v, w);
    results.methods{end+1} = 'GTH';
    results.errors(end+1) = error;
    results.times(end+1) = time;
    results.solver_type{end+1} = 'GTH (uses v,w)';
end

if method_flags.lu
    fprintf('Running Standard LU (does NOT use v,w)...\n');
    [time, error] = mmatrix.run_standard_lu(A, b, x_exact_lu);
    results.methods{end+1} = 'Standard LU';
    results.errors(end+1) = error;
    results.times(end+1) = time;
    results.solver_type{end+1} = 'LU (no v,w)';
end

if method_flags.dmmtf2
    fprintf('Running DMMTF2 (uses v,w)...\n');
    [time, error] = mmatrix.run_dmmtf2(A, b, x_exact_gth, v, w);
    results.methods{end+1} = 'DMMTF2';
    results.errors(end+1) = error;
    results.times(end+1) = time;
    results.solver_type{end+1} = 'GTH (uses v,w)';
end

if method_flags.dmmtrf
    fprintf('Running DMMTRF (nb=%d, uses v,w)...\n', params.nb);
    [time, error] = mmatrix.run_dmmtrf(A, b, x_exact_gth, v, w, params.nb);
    results.methods{end+1} = 'DMMTRF';
    results.errors(end+1) = error;
    results.times(end+1) = time;
    results.solver_type{end+1} = 'GTH (uses v,w)';
end

if method_flags.dmmtrff
    fprintf('Running DMMTRFF (nb=%d, uses v,w)...\n', params.nb);
    [time, error] = mmatrix.run_dmmtrff(A, b, x_exact_gth, v, w, params.nb);
    results.methods{end+1} = 'DMMTRFF';
    results.errors(end+1) = error;
    results.times(end+1) = time;
    results.solver_type{end+1} = 'GTH (uses v,w)';
end

% Display results table
mmatrix.display_single_results_table(results, cond(A));
end

%--------------------------------------------------------------------------
function run_predefined_test(method_flags, params)
% Run dimension study with FEM problem and display tables


% Get parameters
params.sigma = input('Enter sigma value (default: 0.001): ');
if isempty(params.sigma)
    params.sigma = 1e-3;
end

ns_input = input('Enter matrix dimensions as array (default: [128 256 512 1024]): ', 's');
if isempty(ns_input)
    params.ns = [128 256 512 1024];
else
    params.ns = eval(ns_input);
end

% Run experiments
results = mmatrix.run_dimension_experiments(params, method_flags);

% Display results tables
mmatrix.display_dimension_tables(results, params);
end
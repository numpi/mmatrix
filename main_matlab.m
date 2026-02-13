function compare_factorization_methods()
% COMPARE_FACTORIZATION_METHODS Compare different factorization methods for FEM matrices
%
%   Methods compared:
%   - GTH: pure MATLAB's GTH algorithm
%   - Standard LU: MATLAB's built-in LU factorization (This uses dgetrf.f)
%   - dmmtf2: Unblocked version of the GTH algorithm
%   - dmmtrf: Blocked version of the GTH algorithm with DMMTF2 inside
%   - dmmtrff: Blocked version of the GTH algorithm with DMMTRF2 inside

build_mex() % this file complies fortran codes for the first run.
%% Initialization
clear; clc;
fprintf('The condition number varies with σ = 10^{-k} for k = 1 to 12.\n\n');

% Problem parameters
params = struct();
params.n = 512;           % Matrix size (n x n)
params.nb = 64;          % Block size for blocked algorithms
params.gamma = 0;        % Boundary condition parameter
params.sigmas = 10.^(-(1:12));  % Range of sigma values (affects conditioning)
params.h = 1 / (params.n + 1);   % Grid spacing
params.v = ones(params.n, 1);    % Vector for modified LU

fprintf('Problem details:\n');
fprintf('  Matrix size: %dx%d\n', params.n, params.n);
fprintf('  Block size: %d\n', params.nb);
fprintf('  Sigma range: 10^-1 to 10^-%d\n', length(params.sigmas));
fprintf('  Grid spacing h = %.4f\n\n', params.h);

%% Run experiments
results = run_experiments(params);

%% Create figure for outputs
create_visualization(results, params);

end

function results = run_experiments(params)
% RUN_EXPERIMENTS Perform the factorization experiments for all sigma values

% Extract parameters
n = params.n;
nb = params.nb;
gamma = params.gamma;
sigmas = params.sigmas;
h = params.h;
v = params.v;

% Initialize result arrays
num_tests = length(sigmas);
c = 1:num_tests;

results = struct();
results.c = c;
results.conds = zeros(1, num_tests);
results.norms_gth = zeros(1, num_tests);
results.norms_lu = zeros(1, num_tests);
results.norms_dmmtf2 = zeros(1, num_tests);
results.norms_dmmtrf = zeros(1, num_tests);
results.norms_dmmtrff = zeros(1, num_tests);
results.tgth = zeros(1, num_tests);
results.tlu = zeros(1, num_tests);
results.tdmmtf2 = zeros(1, num_tests);
results.tdmmtrf = zeros(1, num_tests);
results.tdmmtrff = zeros(1, num_tests);


for k = c

    sigma = sigmas(k);

    % Create matrix and right-hand side
    [A, b, w] = create_fem_problem(n, h, sigma, gamma);

    % Compute exact solution using high precision
    x_exact = exact_solution(A, v, w, b);
    results.conds(k) = cond(A);

    % Run each factorization method
    [results.tgth(k), results.norms_gth(k), ~] = run_gthlu(A, b, x_exact, v, w);
    [results.tlu(k), results.norms_lu(k), ~] = run_standard_lu(A, b, x_exact);
    [results.tdmmtf2(k), results.norms_dmmtf2(k), info2] = run_dmmtf2(A, b, x_exact, v, w);
    [results.tdmmtrf(k), results.norms_dmmtrf(k), info3] = run_dmmtrf(A, b, x_exact, v, h, sigma, nb);
    [results.tdmmtrff(k), results.norms_dmmtrff(k), info4] = run_dmmtrff(A, b, x_exact, v, h, sigma, nb);

    % Check for warnings
    if info2 ~= 0
        fprintf('\n  Warning: dmmtf2 failed with info=%d at sigma index %d\n', info2, k);
    end
    if info3 ~= 0
        fprintf('\n  Warning: dmmtrf failed with info=%d at sigma index %d\n', info3, k);
    end
    if info4 ~= 0
        fprintf('\n  Warning: dmmtrff failed with info=%d at sigma index %d\n', info4, k);
    end

end

end

function [A, b, w] = create_fem_problem(n, h, sigma, gamma)
% CREATE_FEM_PROBLEM Create FEM matrix and RHS for 1D Poisson problem
% Create tridiagonal FEM matrix
A = full(gallery('tridiag', n, -1, 2, -1));
A(1,1)     = 1 + h*sigma;
A(end,end) = 1 + h*sigma;

% W 
w = zeros(n,1);
w(1) = h*sigma;
w(end) = h*sigma;

% Right-hand side
grid = h * (1:n)';
b = h * sin(pi * grid);
b(1) = b(1)/2 - h*gamma;
b(end) = b(end)/2 + h*gamma;
end

function x_exact = exact_solution(A, v, w, b)
% EXACT_SOLUTION High-precision exact solution
n = length(A);
Avpa = vpa(A);
Avpa(1:n+1:end) = 0;
f = w - Avpa*v;
Avpa(1:n+1:end) = f ./ v;
x_exact = double(Avpa \ b);
end

function [time, rel_error, info] = run_gthlu(A, b, x_exact, v, w)
% RUN_GTH Run GTH factorization
info = 0;
try
    tic;
    [L, U] = gthlu(A, v, w);
    time = toc;
    x = U \ (L \ b);
    rel_error = norm(x - x_exact) / norm(x_exact);
catch ME
    warning('GTH failed: %s', ME.message);
    time = NaN;
    rel_error = NaN;
    info = -1;
end
end

function [time, rel_error, info] = run_standard_lu(A, b, x_exact)
% RUN_STANDARD_LU Run standard LU factorization
info = 0;
try
    tic;
    [L, U] = lu(A);
    time = toc;
    x = U \ (L \ b);
    rel_error = norm(x - x_exact) / norm(x_exact);
catch ME
    warning('Standard LU failed: %s', ME.message);
    time = NaN;
    rel_error = NaN;
    info = -1;
end
end

function [time, rel_error, info] = run_dmmtf2(A, b, x_exact, v, w)
% RUN_DMMTF2 Run dmmtf2 factorization
n = length(A);
try
    tic;
    [Af, info] = dmmtf2_mex(A, v, w);
    time = toc;
    if info == 0
        x = solve_from_factorization(Af, b, n);
        rel_error = norm(x - x_exact) / norm(x_exact);
    else
        rel_error = NaN;
    end
catch ME
    warning('dmmtf2 failed: %s', ME.message);
    time = NaN;
    rel_error = NaN;
    info = -1;
end
end

function [time, rel_error, info] = run_dmmtrf(A, b, x_exact, v, h, sigma, nb)
% RUN_DMMTRF Run dmmtrf factorization
n = length(A);
w_method = zeros(n,1);
w_method(1) = h*sigma;
w_method(end) = h*sigma;

try
    tic;
    [Af, info] = dmmtrf_mex(A, v, w_method, nb);
    time = toc;
    if info == 0
        x = solve_from_factorization(Af, b, n);
        rel_error = norm(x - x_exact) / norm(x_exact);
    else
        rel_error = NaN;
    end
catch ME
    warning('dmmtrf failed: %s', ME.message);
    time = NaN;
    rel_error = NaN;
    info = -1;
end
end

function [time, rel_error, info] = run_dmmtrff(A, b, x_exact, v, h, sigma, nb)
% RUN_DMMTRFF Run dmmtrff factorization
n = length(A);
w_method = zeros(n,1);
w_method(1) = h*sigma;
w_method(end) = h*sigma;

try
    tic;
    [Af, info] = dmmtrff_mex(A, v, w_method, nb);
    time = toc;
    if info == 0
        x = solve_from_factorization(Af, b, n);
        rel_error = norm(x - x_exact) / norm(x_exact);
    else
        rel_error = NaN;
    end
catch ME
    warning('dmmtrff failed: %s', ME.message);
    time = NaN;
    rel_error = NaN;
    info = -1;
end
end

function x = solve_from_factorization(Af, b, n)
% SOLVE_FROM_FACTORIZATION Solve system using compact factorization
L = tril(Af, -1) + eye(n);
U = triu(Af);
y = L \ b;
x = U \ y;
end


function create_visualization(results, params)

c = results.c;
n = params.n;

% Create figure with your exact specifications
figure('Position', [100 100 1200 500], 'Color', 'white');

% SUBPLOT 1: ERRORS (loglog)
subplot(1,2,1)
loglog(c, results.conds*eps(1), 'b', 'LineWidth', 2.5, 'DisplayName', 'cond(A)\epsilon');
hold on
loglog(c, results.norms_gth,    'dg', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'GTH');
loglog(c, results.norms_lu,     'sr', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'LU');         
loglog(c, results.norms_dmmtf2, '*b', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'dmmtf2');
loglog(c, results.norms_dmmtrf, 'om', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'dmmtrf');
loglog(c, results.norms_dmmtrff,'+k', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'dmmtrff');                  

hold off
xlabel('Sigma index k (\sigma = 10^{-k})', 'FontSize', 12);
ylabel('Relative error', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);

% SUBPLOT 2: CPU TIME
subplot(1,2,2)
semilogy(c, results.tgth,    'dg', 'MarkerSize',12, 'MarkerFaceColor','g', 'LineStyle','none', 'DisplayName','GTH');
hold on
semilogy(c, results.tlu,     'sr', 'MarkerSize',12, 'MarkerFaceColor','r', 'LineStyle','none', 'DisplayName','LU');        
semilogy(c, results.tdmmtf2, '*b', 'MarkerSize',12, 'LineStyle','none', 'DisplayName','dmmtf2');                       
semilogy(c, results.tdmmtrf, 'om', 'MarkerSize',12, 'MarkerFaceColor','m', 'LineStyle','none', 'DisplayName','dmmtrf');  
semilogy(c, results.tdmmtrff,'+k', 'MarkerSize',12, 'LineStyle','none', 'DisplayName','dmmtrff');                        

hold off
xlabel('Sigma index k (\sigma = 10^{-k})', 'FontSize', 12);
ylabel('CPU time (s)', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
 
end
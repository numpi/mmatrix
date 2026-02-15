function main_matlab()
%  MAIN_MATLAB Compare different factorization methods for FEM matrices
%
%   Methods compared:
%   - GTH: pure MATLAB's GTH algorithm
%   - Standard LU: MATLAB's built-in LU factorization (This uses dgetrf.f)
%   - dmmtf2: Unblocked version of the GTH algorithm
%   - dmmtrf: Blocked version of the GTH algorithm with DMMTF2 inside
%   - dmmtrff: Blocked version of the GTH algorithm with DMMTRF2 inside

build_mex() % this file complies fortran codes for the first run.
clc;
fprintf('The condition number varies with σ = 10^{-k} for k = 1 to 12.\n\n');

% Problem parameters
params = struct();
params.n = 256;                  % Matrix size (n x n)
params.nb = 64;                 % Block size for blocked algorithms
params.sigmas = 10.^(-(1:12));  % Range of sigma values (affects conditioning)
params.h = 1 / (params.n + 1);   % Grid spacing
params.v = ones(params.n, 1);    % Vector v

fprintf('Problem details:\n');
fprintf('  Matrix size: %dx%d\n', params.n, params.n);
fprintf('  Block size: %d\n', params.nb);
fprintf('  Sigma range: 10^-1 to 10^-%d\n', length(params.sigmas));
fprintf('  Grid spacing h = %.4f\n\n', params.h);

% Run experiments
results = run_experiments(params);
% Create figure for outputs
create_fig(results, params);
end

function results = run_experiments(params)
% RUN_EXPERIMENTS Perform the factorization experiments for all sigma values

% Extract parameters
n = params.n; nb = params.nb; sigmas = params.sigmas; h = params.h; v = params.v;

% Initialize result arrays
nt = length(sigmas); c = 1:nt;
results = struct();
results.c = c; znt=zeros(1, nt); results.conds = znt;
results.norms_gth = znt; results.norms_lu = znt; results.norms_dmmtf2 = znt;
results.norms_dmmtrf = znt; results.norms_dmmtrff = znt; results.tgth = znt; 
results.tlu = znt; results.tdmmtf2 = znt; results.tdmmtrf = znt; results.tdmmtrff = znt;

for k = c
    sigma = sigmas(k);

    % Create matrix and right-hand side
    [A, b, w] = create_fem_problem(n, h, sigma);
    results.conds(k) = cond(A);
    % Compute exact solution using high precision
    x_exact = exact_solution(A, v, w, b);
    
    wt = zeros(n,1); wt(1) = h*sigma; wt(end) = h*sigma;
    w1 = zeros(n,1); w1(1) = h*sigma; w1(end) = h*sigma;
    % Run each factorization method
    [results.tgth(k), results.norms_gth(k)] = run_gthlu(A, b, x_exact, v, w);
    [results.tlu(k), results.norms_lu(k)] = run_standard_lu(A, b, x_exact);
    [results.tdmmtf2(k), results.norms_dmmtf2(k), info2] = run_dmmtf2(A, b, x_exact, v, w);
    [results.tdmmtrf(k), results.norms_dmmtrf(k), info3] = run_dmmtrf(A, b, x_exact, v, wt, nb);
    [results.tdmmtrff(k), results.norms_dmmtrff(k), info4] = run_dmmtrff(A, b, x_exact, v, w1, nb);

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

function [A, b, w] = create_fem_problem(n, h, sigma)
% Create tridiagonal FEM matrix
A = full(gallery('tridiag', n, -1, 2, -1)); A(1,1) = 1 + h*sigma; A(end,end) = 1 + h*sigma;
% w
w = zeros(n,1); w(1) = h*sigma; w(end) = h*sigma;
% Right-hand side
b = h * sin(pi * h * (1:n)'); b(1) = b(1)/2; b(end) = b(end)/2;
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

function [time, rel_error] = run_gthlu(A, b, x_exact, v, w)
% RUN_GTH Run GTH factorization
tic; [L, U] = gthlu(A, v, w); time = toc;
x = U \ (L \ b); rel_error = norm(x - x_exact) / norm(x_exact);
end


function [time, rel_error] = run_standard_lu(A, b, x_exact)
% RUN_STANDARD_LU Run standard LU (dgetrf)
tic; [L, U] = lu(A); time = toc;
x = U \ (L \ b); rel_error = norm(x - x_exact) / norm(x_exact);
end

function [time, rel_error, info] = run_dmmtf2(A, b, x_exact, v, w)
% RUN_DMMTF2 Run dmmtf2 factorization
tic; [Af, info] = dmmtf2_mex(A, v, w); time = toc;
x = solvelu(Af, b); rel_error = norm(x - x_exact) / norm(x_exact);
end

function [time, rel_error, info] = run_dmmtrf(A, b, x_exact, v, wt, nb)
% RUN_DMMTRF Run dmmtrf factorization
tic; [Af, info] = dmmtrf_mex(A, v, wt, nb); time = toc;
x = solvelu(Af, b); rel_error = norm(x - x_exact) / norm(x_exact);
end

function [time, rel_error, info] = run_dmmtrff(A, b, x_exact, v, wt, nb)
% RUN_DMMTRFF Run dmmtrff factorization
tic; [Af, info] = dmmtrff_mex(A, v, wt, nb); time = toc;
x = solvelu(Af, b); rel_error = norm(x - x_exact) / norm(x_exact);
end

function x = solvelu(A, b)
% SOLVELU Solve system using compact factorization
n = length(A); L = tril(A, -1) + eye(n); U = triu(A); y = L \ b; x = U \ y;
end


function create_fig(results, params)
c = results.c;
% Create figure with your exact specifications
figure('Position', [100 100 1200 500]);

% SUBPLOT 1: ERRORS (loglog)
subplot(1,2,1)
loglog(c, results.conds*eps(1), 'b', 'LineWidth', 2.5, 'DisplayName', 'cond(A)\epsilon');
hold on
loglog(c, results.norms_gth,    'dg', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'GTH');
loglog(c, results.norms_lu,     'sr', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'LU');
loglog(c, results.norms_dmmtf2, '*b', 'MarkerSize', 10, 'DisplayName', 'dmmtf2');
loglog(c, results.norms_dmmtrf, 'om', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'dmmtrf');
loglog(c, results.norms_dmmtrff,'+k', 'MarkerSize', 10, 'DisplayName', 'dmmtrff');
xlabel('Sigma index k (\sigma = 10^{-k})', 'FontSize', 12);
ylabel('Relative error', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);

% SUBPLOT 2: CPU TIME
subplot(1,2,2)
semilogy(c, results.tgth,    'dg', 'MarkerSize',12, 'MarkerFaceColor','g', 'DisplayName','GTH');
hold on
semilogy(c, results.tlu,     'sr', 'MarkerSize',12, 'MarkerFaceColor','r', 'DisplayName','LU');
semilogy(c, results.tdmmtf2, '*b', 'MarkerSize',12, 'DisplayName','dmmtf2');
semilogy(c, results.tdmmtrf, 'om', 'MarkerSize',12, 'MarkerFaceColor','m', 'DisplayName','dmmtrf');
semilogy(c, results.tdmmtrff,'+k', 'MarkerSize',12,'DisplayName','dmmtrff');
xlabel('Sigma index k (\sigma = 10^{-k})', 'FontSize', 12);
ylabel('CPU time (s)', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);

end
% This code compares three matrix inversion methods:
%   1. MATLAB inv: Standard built-in matrix inversion
%   2. inv: Inverse computed from compact LU factors and dmmtri
%   3. inv*: Inverse computed by solving AX = I with dmmtrf_solve

clc;
clear;
close all;

% Parameters
n = 32;                             % Matrix size (n x n)
sigmas = 10.^(-(1:0.5:12));         % Range of sigma values (affects conditioning)
ns = length(sigmas);
conds = zeros(1, ns);

% errors
nrm_inv_matlab = zeros(1, ns); nrm_inv_dmmtri = zeros(1, ns); nrm_inv_solve  = zeros(1, ns);
comp_inv_matlab = zeros(1, ns); comp_inv_dmmtri = zeros(1, ns); comp_inv_solve  = zeros(1, ns);

% times
tinv_matlab = zeros(1, ns); tinv_dmmtri = zeros(1, ns); tinv_solve  = zeros(1, ns);

for k = 1:ns
    sigma = sigmas(k);

    % Generate the FEM Matrix and parameters
    [objM, ~] = gettriplet('fem', size=n, left=0, sigma=sigma, gamma=1, h=1/(n+1));
    A = objM.M;
    u = objM.u;
    v = objM.v;

    % Compute exact solution using High-Precision VPA
    Ainv_exact = double(vpa(objM) \ eye(n));
    conds(k)   = cond(A);

    % MATLAB's inv
    tic;
    Ainv_matlab = inv(A);
    tinv_matlab(k) = toc;

    nrm_inv_matlab(k)  = norm(Ainv_matlab - Ainv_exact, 'fro') / norm(Ainv_exact, 'fro');
    comp_inv_matlab(k) = max(max(abs((Ainv_matlab - Ainv_exact) ./ Ainv_exact)));

    % dmmtri
    tic;
    U = mmex.dmmtrff(A, u, v);  % Compute upper triangular factor
    Ainv_dmmtri = dmmtri_mex(U); % Compute inverse from LU factors
    tinv_dmmtri(k) = toc;

    nrm_inv_dmmtri(k)  = norm(Ainv_dmmtri - Ainv_exact, 'fro') / norm(Ainv_exact, 'fro');
    comp_inv_dmmtri(k) = max(max(abs((Ainv_dmmtri - Ainv_exact) ./ Ainv_exact)));

    % dmmtrf_solve (Solving with Identity Matrix) 
    tic;
    Ainv_solve = mmex.dmmtrf_solve(A, u, v, eye(n));
    tinv_solve(k) = toc;

    nrm_inv_solve(k)  = norm(Ainv_solve - Ainv_exact, 'fro') / norm(Ainv_exact, 'fro');
    comp_inv_solve(k) = max(max(abs((Ainv_solve - Ainv_exact) ./ Ainv_exact)));

    fprintf('Iteration %d/%d: cond(A) = %.2e\n', k, ns, conds(k));
end

figure
loglog(1:ns, conds * eps, '-', 'LineWidth', 2.5, 'DisplayName', 'cond(A) * eps');hold on;
loglog(1:ns, nrm_inv_matlab, 'o',  'Color', [0 0.5 0],   'LineWidth', 2, 'DisplayName', 'MATLAB inv (norm)');
loglog(1:ns, comp_inv_matlab, '--', 'Color', [0 0.5 0],   'LineWidth', 1.5, 'DisplayName', 'MATLAB inv (comp)');
loglog(1:ns, nrm_inv_dmmtri, 's',  'Color', [0.8 0 0.8], 'LineWidth', 2, 'DisplayName', 'inv (norm)');
loglog(1:ns, comp_inv_dmmtri, '--', 'Color', [0.8 0 0.8], 'LineWidth', 1.5, 'DisplayName', 'inv (comp)');
loglog(1:ns, nrm_inv_solve,  '^',  'Color', [0 0.6 0.8], 'LineWidth', 2, 'DisplayName', 'inv*  (norm)');
loglog(1:ns, comp_inv_solve,  '--', 'Color', [0 0.6 0.8], 'LineWidth', 1.5, 'DisplayName', 'inv* (comp)');
xlabel('sigma', 'FontSize', 11);
ylabel('Relative Error', 'FontSize', 11);
lgd = legend('Location', 'best');


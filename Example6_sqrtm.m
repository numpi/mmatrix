% -------------------------------------------------------------------------
% This script evaluates the accuracy of different methods for computing
% the matrix square root as the condition number of the matrix varies.
% It compares MATLAB's sqrtm function versus sqrtm from mmatrix-toolbox.
% -------------------------------------------------------------------------
clc; clear;

% Problem parameters
n = 16;                            % matrix dimension (n x n)
sigma_values = 10.^(-(1:0.2:12));  % range of sigma values (affects conditioning)
ns = length(sigma_values);
conds     = zeros(ns, 1);

% Normwise relative errors
err_nrm  = zeros(ns, 1);
err_nrm1 = zeros(ns, 1); % method with maxit = 20
err_nrm2 = zeros(ns, 1); % method with maxit = 30

% Componentwise errors
err_comp  = zeros(ns, 1);
err_comp1 = zeros(ns, 1);
err_comp2 = zeros(ns, 1);

% times
time_matlab = zeros(ns, 1);
time1       = zeros(ns, 1);
time2       = zeros(ns, 1);

for k = 1:ns
    sigma = sigma_values(k);
    % Generate the mmatrix object
    [M_obj, ~] = gettriplet('fem', size=n, left=0, sigma=sigma, gamma=1, h=1/(n+1));
    conds(k) = cond(double(M_obj));
    % ---------------------------------------------------------------------
    % Exact square root (using VPA)
    sqrt_exact = double(sqrtm(vpa(M_obj)));
    % ---------------------------------------------------------------------
    % MATLAB's built-in sqrtm
    sqrt_matlab = sqrtm(double(M_obj));
    err_nrm(k) = norm(sqrt_matlab - sqrt_exact,'fro') / norm(sqrt_exact,'fro');
    err_comp(k) = max(max(abs((sqrt_matlab - sqrt_exact) ./ sqrt_exact)));
    % ---------------------------------------------------------------------
    % sqrtm from mmatrix-toolbox with maxit = 30
    sqrtm_M = sqrtm(M_obj, maxit=20);
    sqrt1 = sqrtm_M.M;
    err_nrm1(k) = norm(sqrt1 - sqrt_exact,'fro') / norm(sqrt_exact,'fro');
    err_comp1(k) = max(max(abs((sqrt1 - sqrt_exact) ./ sqrt_exact)));
    % ---------------------------------------------------------------------
    % sqrtm from mmatrix-toolbox with maxit = 30
    sqrtmM = sqrtm(M_obj, maxit=30);
    sqrt2 = sqrtmM.M;
    err_nrm2(k) = norm(sqrt2 - sqrt_exact, 'fro') / norm(sqrt_exact, 'fro');
    err_comp2(k) = max(max(abs((sqrt2 - sqrt_exact) ./ sqrt_exact)));
    fprintf('Sigma index %2d/%2d: cond(A) = %10.2e\n', k, ns, conds(k));
end

figure
% Plot condition number
loglog(1:ns, conds * eps, 'c', 'DisplayName', 'cond(M) \cdot \epsilon');
hold on;
loglog(1:ns, err_nrm, '-r', 'DisplayName', 'MATLAB sqrtm (Norm)');
loglog(1:ns, err_comp, '--r', 'DisplayName', 'MATLAB sqrtm (Comp)');
loglog(1:ns, err_nrm1, '-b', 'DisplayName', 'sqrtm: 20 iter (Norm)');
loglog(1:ns, err_comp1, '--b','DisplayName', 'sqrtm: 20 iter (Comp)');
loglog(1:ns, err_nrm2, '-m', 'DisplayName', 'sqrtm: 30 iter (Norm)');
loglog(1:ns, err_comp2, '--m','DisplayName', 'sqrtm: 30 iter (Comp)');
xlabel(' \sigma ');
ylabel('Relative Error');
title('Matrix Square Root: MATLAB vs. mmatrix-toolbox');
legend('Location', 'northwest');

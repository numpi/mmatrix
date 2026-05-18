% This script evaluates the accuracy and efficiency of 'eigmin' function
% from mmatrix-toolbox versus MATLAB's standard 'eig' function for an 
% ill-conditioned M-matrix reported in Example 2. 

% Tables generated:
%   Table 1: Fixes dimension (n=100) and varies delta
%   Table 2: Fixes delta (1e-9) and varies dimension (n)

clear; clc;
% ------------------------------------------------------------
% TABLE 1: Accuracy comparison with varying delta

n = 100;
deltas = [1e-3 1e-6 1e-9 1e-12 1e-18 1e-24 1e-30];

fprintf('\nTABLE 1: Varying delta (n = %d)\n', n);
fprintf('delta    lambda     RelErr(eigmin)  RelErr(min(eig(A)))\n');
fprintf('-------------------------------------------------------\n');

for k = 1:length(deltas)
    delta = deltas(k);
    
    % Construct the P matrix
    P = zeros(n);
    for i = 1:n-1
        P(i, i+1) = 1;
    end
    P(n, 1) = delta;
    
    % Compute exact eigenvalue using VPA
    lambda_exact = double(1 - delta^vpa(1/n));
    
    % Define the triplet representation vectors (u, v)
    u = ones(n, 1);
    v = zeros(n, 1);
    v(n) = 1 - delta;
    
    % Construct the ill-conditioned matrix A
    A = diag((v + P*u) ./ u) - P;
    M = mmatrix(A, u, v);
    
    % eigmin
    [lambda_eigmin, ~, ~, ~] = eigmin(M);
    
    % standard MATLAB eig
    lambda_eig = min(eig(A));
    
    % relative errors
    err_eigmin = abs(lambda_eigmin - lambda_exact) / lambda_exact;
    err_eig   = abs(lambda_eig - lambda_exact) / lambda_exact;
    
    fprintf('%1.1e   %1.2e   %1.2e      %1.2e \n', ...
        delta, lambda_exact, err_eigmin, err_eig);
end

% ------------------------------------------------------------
% TABLE 2: Accuracy and CPU-time comparison with varying dimension (n)
% ------------------------------------------------------------
delta = 1e-9;
ns = [10 20 50 100 500 1000 5000];

fprintf('\n\nTABLE 2: Varying dimension (delta = %1.1e)\n', delta);
fprintf('n      lambda     RelErr(eigmin)  RelErr(min(eig))  Time(eigmin)  Time(min(eig))\n');
fprintf('--------------------------------------------------------------------------------\n');

for k = 1:length(ns)
    n = ns(k);
    
    % Construct the P matrix for current dimension
    P = zeros(n);
    for i = 1:n-1
        P(i, i+1) = 1;
    end
    P(n, 1) = delta;
    
    % Compute exact eigenvalue using VPA
    lambda_exact = double(1 - delta^vpa(1/n));
    
    % Define the triplet representation vectors
    u = ones(n, 1);
    v = zeros(n, 1);
    v(n) = 1 - delta;
    
    % Construct matrix A
    A = diag((v + P*u) ./ u) - P;
    M = mmatrix(A, u, v); % M with right triplet representation
    
    % eigmin
    tic; [lambda_eigmin, ~, ~, ~] = eigmin(M); time_eigmin = toc;
    
    % standard MATLAB eig
    tic; lambda_eig = min(eig(A)); time_eig = toc; 

    % relative errors
    err_eigmin = abs(lambda_eigmin - lambda_exact) / lambda_exact;
    err_eig   = abs(lambda_eig - lambda_exact) / lambda_exact;
    
    fprintf('%-4d   %1.2e   %1.2e      %1.2e      %1.4f        %1.4f \n', ...
        n, lambda_exact, err_eigmin, err_eig, time_eigmin, time_eig);
end

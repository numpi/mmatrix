% Example 4 (Accurate singular value decomposition of an M‑matrix).

% This example illustrates how to compute the singular values of a
% randomly generated ill‑conditioned M‑matrix using the mmatrix-toolbox.

% We compare three ways to obtain the singular values:
%   1) mmatrix svd   – accurate, using the triplet representation
%   2) MATLAB's svd  – standard double‑precision svd
%   3) vpa svd       – high‑precision svd (100 decimal digits)
%
% The results shows that the svd from mmatrix-toolbox remains accurate
% down to the level of sigma1*eps, whereas the standard svd loses all
% accuracy for smaller singular values.

clear;close all; clc;
% Matrix sizes to test
n_values = [20, 30];

figure('Position', [100, 100, 1300, 500]);
for idx = 1:length(n_values)
    n = n_values(idx);
    % ------------------------------------------------------------
    A_off = -rand(n); A_off(1:n+1:end) = 0;

    v = zeros(n,1);
    for i = 1:n
        v(i) = rand() * 10^randi([-40,-20]);
    end
    for i = 1:n
        r = rand(); j = randi([-50,50]); scale = r * 10^j;
        A_off(i,:) = A_off(i,:) * scale; v(i) = v(i) * scale;
    end
    % ------------------------------------------------------------
    % Create u with varying magnitudes
    u = rand(n,1).*exp(5*randn(n,1));
    % ------------------------------------------------------------
    % Construct A_true via vpa
    A_off_vpa = vpa(A_off);
    u_vpa     = vpa(u);
    v_vpa     = vpa(v);
    A_true_vpa = A_off_vpa;
    diag_entries = (v_vpa - A_true_vpa*u_vpa) ./ u_vpa;
    for i = 1:n
        A_true_vpa(i,i) = diag_entries(i);
    end
    A_true = double(A_true_vpa);
    % ------------------------------------------------------------
    fprintf('Matrix condition number:     %.2e\n', cond(A_true));
    fprintf('Initial A*u=v (rel error):   %.2e\n', max(double(abs(A_true_vpa*u_vpa - v_vpa))));
    % ------------------------------------------------------------

    [Uf, Sf, Vf] = mmex.dmmsvd(A_off, u, v);
    % the above fortran code was written based on the below matlab code:

    % SVD via standard GTH + pivoiting
    % [Sf, Uf, Vf, L, D, U_factor, P] = SVD2(A_off, u, v);
    % Factorization residual (to make sure if its correct)
    % A_rec = P' * L * D * U_factor * P;
    % fprintf('Factorization Residual:  %.2e\n', max(max(abs(A_true-A_rec) ./ max(abs(A_true)))));
    %------------------------------------------------------------

    % MATLAB standard SVD
    [Um, Sm, Vm] = svd(A_true, 'econ');
    Sm = sort(diag(Sm), 'descend');

    % High precision svd
    [Uv, Sv, Vv] = svd(A_true_vpa, 'econ');
    Sv = sort(double(diag(Sv)), 'descend');
    subplot(1, length(n_values), idx);
    p1 = semilogy(1:n, Sf, 'rx', 'LineWidth', 2, 'MarkerSize', 10); hold on;
    p2 = semilogy(1:n, Sv, 'bs', 'LineWidth', 2, 'MarkerSize', 12);
    p3 = semilogy(1:n, Sm, 'mo', 'LineWidth', 2, 'MarkerSize', 8);

    sigma1_eps = double(Sv(1)) * eps;
    h_line = yline(sigma1_eps, 'k--', 'LineWidth', 1.5);

    text(n*0.5, sigma1_eps*1e5, ...
        sprintf('\\sigma_1 \\epsilon = %.1e', sigma1_eps), ...
        'FontSize', 12, 'FontWeight', 'bold');

    title(sprintf('Accuracy Comparison (n = %d)', n));
    xlabel('Singular Value Index (i)');
    ylabel('Value (\\sigma_i)');
    grid on;

    legend([p1, p2, p3, h_line], ...
        {'svd (accurate)', ...
        'VPA svd', ...
        'MATLAB svd', ...
        '\\sigma_1 \\epsilon'}, ...
        'Location', 'southwest', 'FontSize', 9);
    set(gca, 'FontSize', 11);
end

% ------------------------------------------------------------
function [sigma, U_final, V_final, L, D, U_factor, P] = SVD2(A_off, u, v)
[L, D, U_factor, P] = gthldu(A_off, u, v);
% P*A*P' = L*D*U_factor  % A = P'*L*D*U_factor*P
X = P' * L;
Y = P' * U_factor';
XD = X * D;
[Q, R, perm_vec] = qr(XD, 0);
P_q = eye(size(R,2));
P_q = P_q(:, perm_vec);
W = R * (P_q' * Y');
[U_W, S, V_W] = svd(W, 'econ');
U_final = Q * U_W;
V_final = V_W;
sigma = sort(diag(S), 'descend');
end
% ------------------------------------------------------------
function [L, D, U_factor, P] = gthldu(M, u, v)
n = length(M);
L = eye(n);
U = M;
P = eye(n);
for k = 1:n
    diags = zeros(n-k+1,1);
    for i = k:n
        s = 0;
        for j = k:n
            if j ~= i
                s = s + U(i,j) * u(j);
            end
        end
        diags(i-k+1) = (v(i) - s) / u(i);
    end
    % Pivot selection
    [~, max_idx] = max(abs(diags));
    p = k + max_idx - 1;
    % pivoting
    if p ~= k
        U([k p], :) = U([p k], :);
        U(:, [k p]) = U(:, [p k]);
        if k > 1
            L([k p],1:k-1) = L([p k],1:k-1);
        end
        u([k p]) = u([p k]);
        v([k p]) = v([p k]);
        P([k p], :) = P([p k], :);
    end
    U(k,k) = diags(max_idx);
    % Elimination
    if k < n
        pivot = U(k,k);
        if pivot ~= 0
            L(k+1:end,k) = U(k+1:end,k) / pivot;
            U(k+1:end,k+1:end) = U(k+1:end,k+1:end) - L(k+1:end,k) * U(k,k+1:end);
            v(k+1:end) = v(k+1:end) - L(k+1:end,k) * v(k);
        end
        U(k+1:end,k) = 0;
    end
end
D = diag(diag(U));
U_factor = eye(n);
for i = 1:n
    if D(i,i) ~= 0
        U_factor(i,i+1:end) =  U(i,i+1:end) / D(i,i);
    end
end
end

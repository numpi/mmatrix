% Example 4 (Accurate singular value decomposition of an M‑matrix).

% This example illustrates how to compute the singular values of a
% randomly generated ill‑conditioned M‑matrix using the mmatrix-toolbox.
% The construction follows the method described in:
%   J. Demmel and P. Koev, "Accurate SVDs of weakly diagonally dominant
%   M‑matrices", Numerische Mathematik, 98:99–104, 2004.
%
% We compare three ways to obtain the singular values:
%   1) mmatrix svd   – accurate, using the triplet representation
%   2) MATLAB's svd  – standard double‑precision svd
%   3) vpa svd       – high‑precision svd (100 decimal digits)
%
% The results shows that the svd from mmatrix-toolbox remains accurate
% down to the level of sigma1*eps, whereas the standard svd loses all 
% accuracy for smaller singular values.

clear; clc; close all;

% Matrix sizes to test
n_values = [20, 40];

figure('Position', [100, 100, 1300, 500]);
for idx = 1:length(n_values)
    n = n_values(idx);
    A_off = -rand(n);               % uniform in [-1,0]
    A_off(1:n+1:end) = 0;           % zero diagonal

    % Row sums s_i = r * 10^k, with r in (0,1) and k in [-40, -20].
    % This makes the matrix ill‑conditioned.
    s = zeros(n,1);
    for i = 1:n
        r = rand();
        k = randi([-40, -20]);
        s(i) = r * 10^k;
    end
    % Additional scaling: multiply each row i and s_i by r * 10^j,
    % where j in [-40, 40]
    for i = 1:n
        r = rand();
        j = randi([-40, 40]);
        scale = r * 10^j;
        A_off(i,:) = A_off(i,:) * scale;
        s(i) = s(i) * scale;
    end
    % Build the true M‑matrix:
    A_true = A_off;
    A_true(1:n+1:end) = s - sum(A_off, 2);
    % -------------------------------------------------------------
    % Compute SVD using the mmatrix-toolbox
    u = ones(n,1); v = s;
    objM = mmatrix(A_off, u, v);  
    [Uf, Sf_mat, Vf] = svd(objM); Sf = diag(Sf_mat); 
    % -------------------------------------------------------------
    % SVD using VPA
    Avpa = vpa(A_true,100);
    Avpa(1:n+1:end) = 0;
    Avpa(1:n+1:end) = (v - Avpa*u) ./ u;
    [Uv, Sv_mat, Vv] = svd(Avpa); Sv = (diag(Sv_mat));
    
    % we can also use mmatrix-toolbox to obtain high-precision SVD:
    % [Uv, Sv_mat, Vv] = svd(vpa(objM,100)); Sv = double(diag(Sv_mat)); 
    % -------------------------------------------------------------
    % MATLAB SVD  
    [Um, Sm_mat, Vm] = svd(A_true);
    Sm = diag(Sm_mat);
    % -------------------------------------------------------------
    % Plot results
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
 

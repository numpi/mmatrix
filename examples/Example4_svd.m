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

clear; close all; clc;

n_values = [20 40];
u_modes = {'Well-conditioned', 'Ill-conditioned'};
digits(110) 

figure('Position', [100, 100, 1400, 900]);
plot_idx = 1;
for row = 1:length(u_modes)
    u_mode = u_modes{row};
    
    for col = 1:length(n_values)
        n = n_values(col);
        
        % ------------------------------------------------------------
        A_off = -rand(n); A_off(1:n+1:end) = 0;
        v = zeros(n,1);
        for i = 1:n
            v(i) = rand() * 10^randi([-40,-20]);
        end
        for i = 1:n
            r = rand(); j = randi([-80,80]); scale = r * 10^j;
            A_off(i,:) = A_off(i,:) * scale; v(i) = v(i) * scale;
        end
        
        % ------------------------------------------------------------
        if strcmp(u_mode, 'Well-conditioned')
            u = rand(n,1).*exp(1*randn(n,1));
        else
            u = rand(n,1).*exp(15*randn(n,1));
        end
        kappa_u = max(u) / min(u);
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
        
        fprintf('\n[%s u, n = %d]\n', u_mode, n);
        fprintf('Matrix condition number:     %.2e\n', cond(A_true));
        fprintf('kappa(diag(u)):              %.2e\n', kappa_u);
        fprintf('A*u-v (rel error):   %.2e\n', max(double(abs(abs(A_true_vpa*u_vpa) - v_vpa))));
        
        % mmatrix-toolbox SVD
        MM = mmatrix(A_off, u, v);
        [Uf, Sf, Vf] = svd(MM); Sf = diag(Sf);
        
        % MATLAB standard SVD
        [Um, Sm, Vm] = svd(A_true, 'econ');
        Sm = sort(diag(Sm), 'descend');
        
        % High precision svd
        [Uv, Sv, Vv] = svd(A_true_vpa, 'econ');
        Sv = sort(double(diag(Sv)), 'descend');


        % componentwise relative 
        rel_err_Uf = max(max( (abs(Uf) - abs(Uv)) ./ max(abs(Uv)) ));
        rel_err_Vf = max(max( (abs(Vf) - abs(Vv)) ./ max(abs(Vv)) ));
        rel_err_Um = max(max( (abs(Um) -  abs(Uv)) ./ max(abs(Uv)) ));
        rel_err_Vm = max(max( (abs(Vm) - abs(Vv)) ./ max(abs(Vv)) ));
        fprintf('  cond(Uf):            %.2e\n', cond(Uf));
        fprintf('  cond(Um):            %.2e\n', cond(Um));
        fprintf('  Error in U:            %.2e\n', rel_err_Uf);
        fprintf('  Error in V:            %.2e\n', rel_err_Vf);
        fprintf('  Error in U:            %.2e\n', rel_err_Um);
        fprintf('  Error in V:            %.2e\n', rel_err_Vm);
        
        % Plotting
        subplot(2, length(n_values), plot_idx);
        p1 = semilogy(1:n, Sf, 'rx', 'LineWidth', 2, 'MarkerSize', 10); hold on;
        p2 = semilogy(1:n, Sv, 'bs', 'LineWidth', 2, 'MarkerSize', 12);
        p3 = semilogy(1:n, Sm, 'mo', 'LineWidth', 2, 'MarkerSize', 8);
        
        sigma1_eps = double(Sv(1)) * eps;
        h_line = yline(sigma1_eps, 'k--', 'LineWidth', 1.5);
        
        text(n*0.4, sigma1_eps*1e3, ...
            sprintf('\\sigma_1 \\epsilon = %.1e', sigma1_eps), ...
            'FontSize', 10, 'FontWeight', 'bold');
            
        title(sprintf('%s u (n = %d)\n\\kappa(diag(u)) = %.2e', u_mode, n, kappa_u));
        xlabel('Singular Value Index (i)');
        ylabel('Value (\sigma_i)');
        grid on;
        
        legend([p1, p2, p3, h_line], ...
            {'svd (accurate)', 'VPA svd', 'MATLAB svd', '\sigma_1 \epsilon'}, ...
            'Location', 'southwest', 'FontSize', 8);
        set(gca, 'FontSize', 10);
        
        plot_idx = plot_idx + 1;
    end
end

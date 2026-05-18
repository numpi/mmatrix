% This code is for Example 4 (svd).
% To compile mex file:
mex -c -largeArrayDims dmmltf.f
mex -c -largeArrayDims dmmsvd.f
mex -largeArrayDims dmmsvd_mex.c dmmltf.obj dmmsvd.obj -lmwlapack -lmwblas

clear;
n_values = [20,40];
figure('Position', [100, 100, 1300, 500]);

for idx = 1:length(n_values)
    n = n_values(idx);
    % generate A
    A_off = -rand(n); A_off(1:n+1:end) = 0;
    s = zeros(n,1);
    for i = 1:n
        r = rand(); k = randi([-40, -20]); s(i) = r * 10^k;
    end
    for i = 1:n
        r = rand(); j = randi([-40, 40]); scale = r * 10^j;
        A_off(i,:) = A_off(i,:) * scale;
        s(i) = s(i) * scale;
    end

    % Build A_true:
    A_true = A_off;
    A_true(1:n+1:end) = s - sum(A_off, 2);

    % Compute SVDs:
    [Uf, Sf, Vf] = dmmsvd_mex(A_off, s);
    [Um, Sm, Vm] = svd(A_true);
    [Uv, Sv, Vv] = svd(vpa(A_true, 200));
    Sv = double(Sv);

    subplot(1, length(n_values), idx);
    p1 = semilogy(1:n, Sf, 'rx', 'LineWidth', 2, 'MarkerSize', 10); hold on;
    p2 = semilogy(1:n, diag(Sv), 'bs', 'LineWidth', 2, 'MarkerSize', 12);
    p3 = semilogy(1:n, diag(Sm), 'mo', 'LineWidth', 2, 'MarkerSize', 8);

    sigma1_eps = Sv(1) * eps;
    h_line = yline(sigma1_eps, 'k--', 'LineWidth', 1.5);

    text(n*0.5, sigma1_eps*1e5, sprintf('\\sigma_1 \\epsilon = %.1e', sigma1_eps), ...
        'FontSize', 12, 'FontWeight', 'bold');

    title(sprintf('Accuracy Comparison (n=%d)', n));
    xlabel('Singular Value Index (i)'); ylabel('Value (\sigma_i)');
    grid on;
    legend([p1, p2, p3, h_line], ...
        {'dmmsvd (MEX)', 'VPA', 'MATLAB svd (Standard)', '\sigma_1\epsilon'}, ...
        'Location', 'southwest', 'FontSize', 9);

    set(gca, 'FontSize', 11);
end

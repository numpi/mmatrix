clear; close all; clc;
digits(300)

n = 10;
rnd=rand(n);r1=rand(n,1);r2=randn(n,1);

% ------------------------------------------------------------
u_scales = [0 1 2 3 4 6 8 9 10 11 12 13 14 15 16 18 20 25];

num_tests = length(u_scales);

kappa_vals = zeros(num_tests,1);
errL_vals  = zeros(num_tests,1);
errU_vals  = zeros(num_tests,1);
errUf_vals  = zeros(num_tests,1);
errUm_vals  = zeros(num_tests,1);
errD_vals  = zeros(num_tests,1);
errVf_vals  = zeros(num_tests,1);
errVm_vals  = zeros(num_tests,1);
errVf  = zeros(num_tests,1);
errVm  = zeros(num_tests,1);
errUf  = zeros(num_tests,1);
errUm  = zeros(num_tests,1);
for t = 1:num_tests
    A_off = -rnd;
    A_off(1:n+1:end) = 0;
    v = zeros(n,1);
    for i = 1:n
        v(i) = rand() * 10^randi([-40,-20]);
    end
    for i = 1:n
        r = rand();
        j = randi([-80,80]);
        scale = r * 10^j;
        A_off(i,:) = A_off(i,:) * scale;
        v(i)       = v(i) * scale;
    end
    % ------------------------------------------------------------
    % Construct u
    alpha = u_scales(t);
    if alpha==0
        u=ones(n,1);
    else
        u  = r1.*exp(alpha*r2);
    end
    kappa_u = max(abs(u)) / min(abs(u));
    % ------------------------------------------------------------
    % Build exact matrix in VPA
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
    % Double precision factorization
    [L,D,U_factor,P] = gthldu(A_off,u,v);
    % High precision
    [Lv,Dv,Uv_ref,Pv] = gthldu(vpa(A_off),vpa(u),vpa(v));

    Lv = double(Lv);
    Dv = double(Dv);
    Uv_ref = double(Uv_ref);
    % ------------------------------------------------------------
    % Componentwise relative errors
    % ------------------------------------------------------------
    errL = abs(L - Lv) ./ (abs(Lv));
    errU = abs(U_factor - Uv_ref) ./ (abs(Uv_ref));
    errD = abs(D - Dv) ./ (abs(Dv));

    errL_vals(t) = max(errL(:));
    errU_vals(t) = max(errU(:));
    errD_vals(t) = max(errD(:));

    kappa_vals(t) = kappa_u;

    % mmatrix-toolbox SVD
    [Uf, Sf, Vf] = mmex.dmmsvd(A_off, u, v); %Or Obj = mmatrix(A_off, u, v);[Uf, Sf, Vf] = svd(Obj); Sf = diag(Sf);

    % MATLAB standard SVD
    [Um, Sm, Vm] = svd(A_true);

    % High precision svd
    [Uv, Sv, Vv] = svd(A_true_vpa);

    Uf=abs(Uf);Uv=abs(Uv);Um=abs(Um);Vf=abs(Vf);Vv=abs(Vv);Vm=abs(Vm);
    % componentwise relative
    e1 = abs(Uf-Uv) ./ max(Uv,eps); errUf_vals(t)=max(e1(:));
    e2 = abs(Vf-Vv) ./ max(Vv,eps);errVf_vals(t)=max(e2(:));
    e3 = abs(Um-Uv) ./ max(Uv,eps);errUm_vals(t) =max(e3(:));
    e4 = abs(Vm-Vv) ./ max(Vv,eps);errVm_vals(t) =max(e4(:));

    errUf(t) = norm(Uf - Uv,'fro')/norm(Uv,'fro');
    errVf(t) = norm(Vf - Vv,'fro')/norm(Vv,'fro');
    errUm(t) = norm(Um -  Uv,'fro')/norm(Uv,'fro');
    errVm(t) = norm(Vm - Vv,'fro')/norm(Vv,'fro');
    fprintf('\nTest %d\n', t);
    fprintf('kappa(diag(u)) = %.2e\n', kappa_u);
    fprintf('cw error L     = %.2e\n', errL_vals(t));
    fprintf('cw error U     = %.2e\n', errU_vals(t));
    fprintf('cw error D     = %.2e\n', errD_vals(t));
    fprintf('cw error Uf     = %.2e\n', errUf_vals(t));
    fprintf('cw error Um     = %.2e\n', errUm_vals(t));
    fprintf('cw error Vf    = %.2e\n', errVf_vals(t));
    fprintf('cw error Vm    = %.2e\n', errVm_vals(t));
    fprintf('nw error Uf     = %.2e\n', errUf(t));
    fprintf('nw error Um     = %.2e\n', errUm(t));
    fprintf('nw error Vf    = %.2e\n', errVf(t));
    fprintf('nw error Vm    = %.2e\n', errVm(t));
end

% ------------------------------------------------------------
% Plot error
% ------------------------------------------------------------
opt = {'LineStyle', 'none', 'MarkerSize', 8, 'LineWidth', 1.5};
loglog(kappa_vals, errL_vals,  'o', opt{:}); hold on;
loglog(kappa_vals, errU_vals,  's', opt{:});
loglog(kappa_vals, errD_vals,  '^', opt{:});
loglog(kappa_vals, errUf_vals, 'd', opt{:});
loglog(kappa_vals, errUm_vals, 'v', opt{:});
loglog(kappa_vals, errUf,      'p', opt{:});
loglog(kappa_vals, errUm,      'h', opt{:});
loglog(kappa_vals, kappa_vals*eps, 'm--','LineWidth',2);
xlabel('\kappa(diag(u))');
ylabel('Maximum componentwise relative error');
legend('L, ec','U, ec','D, ec','U fortran SVD, ec', 'U matlab SVD, ec','U fortran SVD, en','U matlab SVD, en','\kappa(diag(u)) \cdot eps','Location','northwest');
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

function main_matlab()
%  Solves a linear system with different factorization methods
%
%   Methods compared:
%   - GTH: pure MATLAB's GTH algorithm
%   - Standard LU: MATLAB's built-in LU factorization (This uses dgetrf.f)
%   - Left  triplet methods (dmltf2, dmltrf2, dmltrf, dmltrff)
%   - Right triplet methods (dmmtf2, dmmtrf2, dmmtrf, dmmtrff)
%   - dmmtf2, dmltf2: Unblocked versions of the GTH algorithm
%   - dmltrf, dmmtrf: Blocked versions using dmltf2 and dmmtf2 inside
%   - dmltrff, dmmtrff: Blocked versions with dmltrf2 and dmmtrf2 inside

% build_mex() % this complies fortran codes

clc;
fprintf('The condition number varies with σ = 10^{-k} for k = 1 to 12.\n\n');
 
% Problem parameters
params.n = 128;                   % Matrix size (n x n)
params.sigmas = 10.^(-(1:12));   % Range of sigma values (affects conditioning)
params.h = 1 / (params.n + 1);   % Grid spacing
params.v = ones(params.n, 1);    % Vector v

fprintf('Problem details:\n');
fprintf('  Matrix size: %dx%d\n', params.n, params.n);
fprintf('  Sigma range: 10^-1 to 10^-%d\n', length(params.sigmas));
fprintf('  Grid spacing h = %.4f\n\n', params.h);

results = run_experiments(params);
create_fig(results);
end

function r = run_experiments(params)

n = params.n; sigmas = params.sigmas; h = params.h; v = params.v;

nt = length(sigmas);
z = zeros(1,nt);

r.c = 1:nt; r.conds = z;

% Errors
r.nrm_gth    = z; r.nrm_lu      = z;
r.nrm_dmmtf2 = z; r.nrm_dmmtrf2 = z; r.nrm_dmmtrf = z; r.nrm_dmmtrff = z;
r.nrm_dmltf2 = z; r.nrm_dmltrf2 = z; r.nrm_dmltrf = z; r.nrm_dmltrff = z;

% Times
r.tgth    = z; r.tlu      = z;
r.tdmmtf2 = z; r.tdmmtrf2 = z; r.tdmmtrf = z; r.tdmmtrff = z;
r.tdmltf2 = z; r.tdmltrf2 = z; r.tdmltrf = z; r.tdmltrff = z;

for k = 1:nt

    sigma = sigmas(k);

    [A,b,w] = create_fem_problem(n,h,sigma);
    r.conds(k) = cond(A);

    x_exact = exact_solution(A,v,w,b);

    % GTH with right triplet
    tic; [L,U] = gthlu(A,v,w); x = U\(L\b); r.tgth(k) = toc;
    r.nrm_gth(k) = norm(x-x_exact)/norm(x_exact);

    % LU
    tic; [L,U] = lu(A); x = U\(L\b); r.tlu(k) = toc;
    r.nrm_lu(k) = norm(x-x_exact)/norm(x_exact);

    % RIGHT TRIPLET--------------------------------------------------------
    tic; x = dmmtf2_mexsolve(A,v,w,b); r.tdmmtf2(k) = toc;
    r.nrm_dmmtf2(k) = norm(x-x_exact)/norm(x_exact);

    tic;  x = dmmtrf2_mexsolve(A,v,w,b); r.tdmmtrf2(k) = toc;
    r.nrm_dmmtrf2(k) = norm(x-x_exact)/norm(x_exact);

    tic;  x = dmmtrf_mexsolve(A,v,w,b); r.tdmmtrf(k) = toc;
    r.nrm_dmmtrf(k) = norm(x-x_exact)/norm(x_exact);

    tic; x = dmmtrff_mexsolve(A,v,w,b); r.tdmmtrff(k) = toc;
    r.nrm_dmmtrff(k) = norm(x-x_exact)/norm(x_exact);

    %LEFT TRIPLET ---------------------------------------------------------
    tic; x = dmltf2_mexsolve(A,v,w,b); r.tdmltf2(k) = toc;
    r.nrm_dmltf2(k) = norm(x-x_exact)/norm(x_exact);

    tic; x = dmltrf2_mexsolve(A,v,w,b); r.tdmltrf2(k) = toc;
    r.nrm_dmltrf2(k) = norm(x-x_exact)/norm(x_exact);

    tic; x = dmltrf_mexsolve(A,v,w,b); r.tdmltrf(k) = toc;
    r.nrm_dmltrf(k) = norm(x-x_exact)/norm(x_exact);

    tic; x = dmltrff_mexsolve(A,v,w,b); r.tdmltrff(k) = toc;
    r.nrm_dmltrff(k) = norm(x-x_exact)/norm(x_exact);
end
end

% PROBLEM -----------------------------------------------------------------
function [A,b,w] = create_fem_problem(n,h,sigma)
A = full(gallery('tridiag',n,-1,2,-1));
A(1,1) = 1 + h*sigma;
A(end,end) = 1 + h*sigma;

w = zeros(n,1);
w(1) = h*sigma;
w(end) = h*sigma;

b = h*sin(pi*h*(1:n)');
b(1) = b(1)/2;
b(end) = b(end)/2;
end

function x = exact_solution(A,v,w,b)
n = length(A);
Avpa = vpa(A);
Avpa(1:n+1:end) = 0;
f = w - Avpa*v;
Avpa(1:n+1:end) = f ./ v;
x = double(Avpa\b);
end


function create_fig(results)
c = results.c;
figure('Position',[100 100 1300 500])

subplot(1,2,1)
loglog(c,results.conds*eps,'c','LineWidth', 2.5,'LineWidth',3,'DisplayName','cond(A)\epsilon')
hold on
loglog(c,results.nrm_lu,'.r','LineWidth',1,'MarkerSize',11,'DisplayName','LU(dgetrf)')
loglog(c,results.nrm_gth,'sr','LineWidth',1,'MarkerSize',7,'DisplayName','GTH(right)')
loglog(c,results.nrm_dmmtf2,'ob','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtf2')
loglog(c,results.nrm_dmmtrf2,'db','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtrf2')
loglog(c,results.nrm_dmmtrf,'+b','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtrf')
loglog(c,results.nrm_dmmtrff,'*b','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtrff')
loglog(c,results.nrm_dmltf2,'om','LineWidth',1,'MarkerSize',9,'DisplayName','dmltf2')
loglog(c,results.nrm_dmltrf2,'dm','LineWidth',1,'MarkerSize',9,'DisplayName','dmltrf2')
loglog(c,results.nrm_dmltrf,'+m','LineWidth',1,'MarkerSize',9,'DisplayName','dmltrf')
loglog(c,results.nrm_dmltrff,'*m','LineWidth',1,'MarkerSize',9,'DisplayName','dmltrff')
xlabel('Sigma','FontSize', 12);
ylabel('Relative error','FontSize', 12);
legend('Location','best','FontSize', 12);
subplot(1,2,2)
semilogy(c,results.tlu,'.r','LineWidth',1,'MarkerSize',11,'DisplayName','LU(dgetrf)')
hold on
semilogy(c,results.tgth,'sr','LineWidth',1,'MarkerSize',7,'DisplayName','GTH(right)')
semilogy(c,results.tdmmtf2,'ob','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtf2')
semilogy(c,results.tdmmtrf2,'db','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtrf2')
semilogy(c,results.tdmmtrf,'+b','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtrf')
semilogy(c,results.tdmmtrff,'*b','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtrff')
semilogy(c,results.tdmltf2,'om','LineWidth',1,'MarkerSize',9,'DisplayName','dmltf2')
semilogy(c,results.tdmltrf2,'dm','LineWidth',1,'MarkerSize',9,'DisplayName','dmltrf2')
semilogy(c,results.tdmltrf,'+m','LineWidth',1,'MarkerSize',9,'DisplayName','dmltrf')
semilogy(c,results.tdmltrff,'*m','LineWidth',1,'MarkerSize',9,'DisplayName','dmltrff')
xlabel('Sigma','FontSize', 12);
ylabel('CPU time (s)','FontSize', 12);
legend('Location','best','FontSize', 12);
end

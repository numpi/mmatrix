function main_matlab()
%  MAIN_MATLAB Compare different factorization methods for FEM matrices
%
%   Methods compared:
%   - GTH: pure MATLAB's GTH algorithm
%   - Standard LU: MATLAB's built-in LU factorization (This uses dgetrf.f)
%   - Left  triplet methods (dmltf2, dmltrf, dmltrff)
%   - Right triplet methods (dmmtf2, dmmtrf, dmmtrff)
%   - dmmtf2,dmltf2: Unblocked version of the GTH algorithm
%   - dmltrf,dmmtrf: Blocked version of the GTH algorithm with DMMTF2 inside
%   - dmltrff,dmmtrff: Blocked version of the GTH algorithm with DMMTRF2 inside

build_mex() % this complies fortran codes 
clc;
fprintf('The condition number varies with σ = 10^{-k} for k = 1 to 12.\n\n');

% Problem parameters
params.n = 256;                  % Matrix size (n x n)
params.nb = 64;                  % Block size for blocked algorithms
params.sigmas = 10.^(-(1:12));   % Range of sigma values (affects conditioning)
params.h = 1 / (params.n + 1);   % Grid spacing
params.v = ones(params.n, 1);    % Vector v

fprintf('Problem details:\n');
fprintf('  Matrix size: %dx%d\n', params.n, params.n);
fprintf('  Block size: %d\n', params.nb);
fprintf('  Sigma range: 10^-1 to 10^-%d\n', length(params.sigmas));
fprintf('  Grid spacing h = %.4f\n\n', params.h);

% Run experiments
results = run_experiments(params);
create_fig(results);

 end

function results = run_experiments(params)
% RUN_EXPERIMENTS Perform the factorization experiments for all sigma values

% Extract parameters
n = params.n; nb = params.nb; sigmas = params.sigmas; h = params.h; v = params.v;


nt = length(sigmas);z = zeros(1,nt);

results.c = 1:nt;results.conds = z;

% Errors
results.norms_gth     = z;results.norms_lu      = z;
results.norms_dmmtf2  = z;results.norms_dmmtrf  = z;results.norms_dmmtrff  = z;
results.norms_dmltrf2 = z;results.norms_dmltrf  = z;results.norms_dmltrff  = z;

% Times
results.tgth     = z;results.tlu     = z;
results.tdmmtf2  = z;results.tdmmtrf = z;results.tdmmtrff = z;
results.tdmltrf2 = z;results.tdmltrf = z;results.tdmltrff  = z;

for k = 1:nt

    sigma = sigmas(k);
    % Create matrix and right-hand side
    [A,b,w] = create_fem_problem(n,h,sigma);
    results.conds(k) = cond(A);
    % Compute exact solution using high precision
    x_exact = exact_solution(A,v,w,b);

    % RIGHT TRIPLET
    [results.tdmmtf2(k), results.norms_dmmtf2(k), ~] = run_dmmtf2(A,b,x_exact,v,w);
    [results.tdmmtrf(k), results.norms_dmmtrf(k), ~] = run_dmmtrf(A,b,x_exact,v,w,nb);
    [results.tdmmtrff(k), results.norms_dmmtrff(k), ~] = run_dmmtrff(A,b,x_exact,v,w,nb);

    % LEFT TRIPLET
    [results.tdmltf2(k), results.norms_dmltf2(k), ~] = run_dmltf2(A,b,x_exact,v,w);
    [results.tdmltrf(k), results.norms_dmltrf(k), ~] = run_dmltrf(A,b,x_exact,v,w,nb);
    [results.tdmltrff(k), results.norms_dmltrff(k), ~] =  run_dmltrff_left(A,b,x_exact,v,w,nb);

    [results.tgth(k), results.norms_gth(k), ~] = run_gth(A,b,x_exact,v,w);
    [results.tlu(k), results.norms_lu(k), ~] = run_lu(A,b,x_exact);

end
end

% PROBLEM 
function [A,b,w] = create_fem_problem(n,h,sigma)
% Create tridiagonal FEM matrix
A = full(gallery('tridiag',n,-1,2,-1));
A(1,1) = 1 + h*sigma;
A(end,end) = 1 + h*sigma;
% this w also works for the left triplet
w = zeros(n,1);w(1) = h*sigma;w(end) = h*sigma; 
% Right-hand side
b = h*sin(pi*h*(1:n)');
b(1) = b(1)/2;
b(end) = b(end)/2;

end


function x = exact_solution(A,v,w,b)
% High precision reference solution
n = length(A);
Avpa = vpa(A);
Avpa(1:n+1:end) = 0;
f = w - Avpa*v;
Avpa(1:n+1:end) = f ./ v;
x = double(Avpa\b);
end

function [time, rel_error, info] = run_gth(A,b,x_exact,v,w)
tic; [L,U] = gthlu(A,v,w); time = toc;x = U\(L\b);
rel_error = norm(x-x_exact)/norm(x_exact);info = 0;
end

function [time, rel_error, info] = run_lu(A,b,x_exact)
tic; [L,U] = lu(A); time = toc;
x = U\(L\b);rel_error = norm(x-x_exact)/norm(x_exact);info = 0;
end


%RIGHT TRIPLET-----------------------------------------------
function [time, rel_error, info] = run_dmmtf2(A,b,x_exact,v,w)
% RUN_DMMTF2
tic; [Af,info] = dmmtf2_mex(A,v,w); time = toc;
x = solvelu(Af,b);
rel_error = norm(x-x_exact)/norm(x_exact);
end

function [time, rel_error, info] = run_dmmtrf(A,b,x_exact,v,w,nb)
% RUN_DMMTRF
tic; [Af,info] = dmmtrf_mex(A,v,w,nb); time = toc;
x = solvelu(Af,b);
rel_error = norm(x-x_exact)/norm(x_exact);
end

function [time, rel_error, info] = run_dmmtrff(A,b,x_exact,v,w,nb)
% RUN_DMMTRFF
tic; [Af,info] = dmmtrff_mex(A,v,w,nb); time = toc;
x = solvelu(Af,b);rel_error = norm(x-x_exact)/norm(x_exact);
end

% LEFT TRIPLET-------------------------------------------------
function [time, rel_error, info] = run_dmltf2(A,b,x_exact,v,w)
% RUN_DMLTRF2
tic; [Af,info] = dmltf2_mex(A,v,w); time = toc;
x = solvelu(Af,b);rel_error = norm(x-x_exact)/norm(x_exact);
end

function [time, rel_error, info] = run_dmltrf(A,b,x_exact,v,w,nb)
% RUN_DMLTRF
tic; [Af,info] = dmltrf_mex(A,v,w,nb); time = toc;
x = solvelu(Af,b);rel_error = norm(x-x_exact)/norm(x_exact);
end

function [time, rel_error, info] = run_dmltrff_left(A,b,x_exact,v,w,nb)
% RUN_DMLTRFF (Left)
tic; [Af,info] = dmltrff_mex(A,v,w,nb); time = toc;
x = solvelu(Af,b);rel_error = norm(x-x_exact)/norm(x_exact);
end


function x = solvelu(A,b)
% Solve compact LU 
n = length(A);L = tril(A,-1) + eye(n);U = triu(A);y = L\b;x = U\y;
end


function create_fig(results)
c = results.c;
figure('Position',[100 100 1300 500])

subplot(1,2,1)
loglog(c,results.conds*eps,'c','LineWidth', 2.5,'LineWidth',3,'DisplayName','cond(A)\epsilon')
hold on
loglog(c,results.norms_lu,'.r','LineWidth',1,'MarkerSize',11,'DisplayName','LU(dgetrf)')
loglog(c,results.norms_gth,'sr','LineWidth',1,'MarkerSize',7,'DisplayName','GTH(right)')
loglog(c,results.norms_dmmtf2,'ob','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtf2')
loglog(c,results.norms_dmmtrf,'+b','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtrf')
loglog(c,results.norms_dmmtrff,'*b','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtrff')
loglog(c,results.norms_dmltf2,'om','LineWidth',1,'MarkerSize',9,'DisplayName','dmltf2')
loglog(c,results.norms_dmltrf,'+m','LineWidth',1,'MarkerSize',9,'DisplayName','dmltrf')
loglog(c,results.norms_dmltrff,'*m','LineWidth',1,'MarkerSize',9,'DisplayName','dmltrff')
xlabel('Sigma','FontSize', 12);
ylabel('Relative error','FontSize', 12);
legend('Location','best','FontSize', 12);
subplot(1,2,2)
semilogy(c,results.tlu,'.r','LineWidth',1,'MarkerSize',11,'DisplayName','LU(dgetrf)')
hold on
semilogy(c,results.tgth,'sr','LineWidth',1,'MarkerSize',7,'DisplayName','GTH(right)')
semilogy(c,results.tdmmtf2,'ob','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtf2')
semilogy(c,results.tdmmtrf,'+b','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtrf')
semilogy(c,results.tdmmtrff,'*b','LineWidth',1,'MarkerSize',9,'DisplayName','dmmtrff')
semilogy(c,results.tdmltf2,'om','LineWidth',1,'MarkerSize',9,'DisplayName','dmltf2')
semilogy(c,results.tdmltrf,'+m','LineWidth',1,'MarkerSize',9,'DisplayName','dmltrf')
semilogy(c,results.tdmltrff,'*m','LineWidth',1,'MarkerSize',9,'DisplayName','dmltrff')
xlabel('Sigma','FontSize', 12);
ylabel('CPU time (s)','FontSize', 12);
legend('Location','best','FontSize', 12);
end
function main_matlab_dimensions()
% Compare factorization methods for ONE sigma and varying dimensions

fprintf('Fixed sigma, varying matrix dimension n.\n\n');

%  ------------------- PARAMETERS -------------------
params = struct();
params.sigma = 1e-9;                 
params.ns = [128 256 512 1024 ];     % Different matrix sizes
params.nb = 64;                      % Block size

% ------------------- RUN EXPERIMENTS --------
results = run_experiments(params);
% ------------------- PLOT -------------------
create_fig(results, params);
end
%----------------------------------------------
function results = run_experiments(params)
sigma = params.sigma;ns = params.ns;nb = params.nb;

nt = length(ns);
results.ns = ns;z = zeros(1, nt);
results.conds = z;results.norms_gth = z;results.norms_lu = z;results.norms_dmmtf2 = z;
results.norms_dmmtrf = z;results.norms_dmmtrff = z;
results.tgth = z;results.tlu = z;results.tdmmtf2 = z;results.tdmmtrf = z;results.tdmmtrff = z;

for k = 1:nt
    n = ns(k);
    h = 1/(n+1);
    v = ones(n,1);

    fprintf('Running n = %d\n', n);

    % Create problem
    [A, b, w] = create_fem_problem(n, h, sigma);

    results.conds(k) = cond(A);

    % Exact solution
    x_exact = exact_solution(A, v, w, b);

    % Run methods
    [results.tgth(k), results.norms_gth(k)] = run_gthlu(A, b, x_exact, v, w);
    [results.tlu(k), results.norms_lu(k)] =  run_standard_lu(A, b, x_exact);
    [results.tdmmtrf(k), results.norms_dmmtrf(k), info3] =  run_dmmtrf(A, b, x_exact, v, w, nb);
    [results.tdmmtrff(k), results.norms_dmmtrff(k), info4] =  run_dmmtrff(A, b, x_exact, v, w, nb);
    [results.tdmmtf2(k), results.norms_dmmtf2(k), info2] =  run_dmmtf2(A, b, x_exact, v, w);

    % Warnings
    if info2 ~= 0
        fprintf('  Warning: dmmtf2 failed (info=%d)\n', info2);
    end
    if info3 ~= 0
        fprintf('  Warning: dmmtrf failed (info=%d)\n', info3);
    end
    if info4 ~= 0
        fprintf('  Warning: dmmtrff failed (info=%d)\n', info4);
    end
end
end
%------------------------------------------------------------------- 
function [A, b, w] = create_fem_problem(n, h, sigma)
A = full(gallery('tridiag', n, -1, 2, -1));A(1,1) = 1 + h*sigma;A(end,end) = 1 + h*sigma;
w = zeros(n,1);w(1) = h*sigma;w(end) = h*sigma;
b = h * sin(pi * h * (1:n)');b(1) = b(1)/2;b(end) = b(end)/2;
end
%------------------------------------------------------------------- 

function x_exact = exact_solution(A, v, w, b)
n = length(A);
Avpa = vpa(A);
Avpa(1:n+1:end) = 0;
f = w - Avpa*v;
Avpa(1:n+1:end) = f ./ v;
x_exact = double(Avpa \ b);
end

%------------------------------------------------------------------- 
function [time, rel_error] = run_gthlu(A, b, x_exact, v, w)
tic;[L, U] = lu_matlab(A, v, w);time = toc;
x = U \ (L \ b);rel_error = norm(x - x_exact) / norm(x_exact);
end

function [time, rel_error] = run_standard_lu(A, b, x_exact)
tic;[L, U] = lu(A);time = toc;
x = U \ (L \ b);rel_error = norm(x - x_exact) / norm(x_exact);
end

function [time, rel_error, info] = run_dmmtf2(A, b, x_exact, v, w)
tic;[Af, info] = mmex.dmmtf2(A, v, w);time = toc;
x = solvelu(Af, b);rel_error = norm(x - x_exact) / norm(x_exact);
end

function [time, rel_error, info] = run_dmmtrf(A, b, x_exact, v, w, nb)
tic;[Af, info] = mmex.dmmtrf(A, v, w, nb);time = toc;
x = solvelu(Af, b);rel_error = norm(x - x_exact) / norm(x_exact);
end

function [time, rel_error, info] = run_dmmtrff(A, b, x_exact, v, w, nb)
tic;[Af, info] = mmex.dmmtrff(A, v, w, nb);time = toc;
x = solvelu(Af, b);rel_error = norm(x - x_exact) / norm(x_exact);
end

function x = solvelu(A, b)
n = length(A);L = tril(A,-1) + eye(n);U = triu(A);y = L \ b;x = U \ y;
end
%------------------------------------------------------------------ 

function create_fig(results, params)
ns = results.ns;
figure('Position', [100 100 1200 500]);
% -------- ERRORS --------
subplot(1,2,1)
loglog(ns, results.conds .* eps, 'k', 'LineWidth',2, ...
    'DisplayName','cond(A)\epsilon'); hold on
loglog(ns, results.norms_gth,    'dg','MarkerSize', 10,'MarkerFaceColor','g','DisplayName','GTH');
loglog(ns, results.norms_lu,     'sr','MarkerSize', 10,'MarkerFaceColor','r','DisplayName','LU');
loglog(ns, results.norms_dmmtf2, '*b','MarkerSize', 10,'DisplayName','dmmtf2');
loglog(ns, results.norms_dmmtrf, 'om','MarkerSize', 10,'MarkerFaceColor','m','DisplayName','dmmtrf');
loglog(ns, results.norms_dmmtrff,'+k','MarkerSize', 10,'DisplayName','dmmtrff');
xlabel('Matrix dimension n');
ylabel('Relative error');
legend('Location','best');
title(sprintf('Relative Errors (\\sigma = %.1e)', params.sigma));
% -------- TIME --------
subplot(1,2,2)
loglog(ns, results.tgth,    'dg', 'MarkerSize',12,'MarkerFaceColor','g','DisplayName','GTH'); hold on
loglog(ns, results.tlu,     'sr', 'MarkerSize',12,'MarkerFaceColor','r','DisplayName','LU');
loglog(ns, results.tdmmtf2, '*b', 'MarkerSize',12,'DisplayName','dmmtf2');
loglog(ns, results.tdmmtrf, 'om', 'MarkerSize',12,'MarkerFaceColor','m','DisplayName','dmmtrf');
loglog(ns, results.tdmmtrff,'+k', 'MarkerSize',12,'DisplayName','dmmtrff');
xlabel('Matrix dimension n');
ylabel('CPU time (s)');
legend('Location','best');
title(sprintf('CPU Time (\\sigma = %.1e)', params.sigma));

end

% QBD: comparison of LRQBD and accLRQBD
%
% Example 5.1 in:
% % G. Gu, L. Wang and R.-C. Li,
% "Highly Accurate Latouche-Ramaswami Logarithmic Reduction
% Algorithm for Quasi-Birth-and-Death Process,"
% Journal of Mathematical Study, 55 (2022), pp. 180-194.

clear; close all;
digits(100);

n = 24; maxit = 40;

% Construct the coefficients in 100-digit arithmetic
r = vpa('1/300');alpha = vpa('18.244');beta = vpa('65536');d = vpa('0.280');

P0 = vpa(zeros(n));P1 = vpa(zeros(n));P2 = vpa(zeros(n));
for i = 1:n
    P0(i,i) = 8*(25-i);
    P2(i,i) = 192*d;
    if i < n
        P1(i,i+1) = alpha*r*(beta-i+1)/beta;
    end
    if i > 1
        P1(i,i-1) = (i-1)*r;
    end
    P1(i,i) = -sum(P1(i,:)) - P0(i,i) - P2(i,i);
end

% Here A1 = 0 and (A0 + A2)*ones(n,1) = ones(n,1).
A0hp = (-P1) \ P0;A2hp = (-P1) \ P2;

% Compute the high-precision solution
Gref = mmatrix.qbd(A0hp, A2hp,maxit=80, accurate=false);

% Run both algorithms in double precision
A0 = double(A0hp);A2 = double(A2hp);
[~, historyLR] = mmatrix.qbd(A0, A2,maxit=maxit, accurate=false);
[~, historyGTH] = mmatrix.qbd(A0, A2, maxit=maxit, accurate=true);
% Compute errors
errorLR  = zeros(maxit+1,2);errorGTH = zeros(maxit+1,2);
for j = 1:maxit+1
    Glr  = vpa(sym(historyLR{j},'f'),100); Ggth = vpa(sym(historyGTH{j},'f'),100);
    Dlr  = abs(Glr-Gref); Dgth = abs(Ggth-Gref);
    errorLR(j,1)  = double(max(Dlr(:)./abs(Gref(:))));
    errorLR(j,2)  = double(norm(Dlr,'fro')/norm(Gref,'fro'));
    errorGTH(j,1) = double(max(Dgth(:)./abs(Gref(:))));
    errorGTH(j,2) = double(norm(Dgth,'fro')/norm(Gref,'fro'));
end
% Plot componentwise and normwise errors
figure;
k = 0:maxit;
semilogy(k,errorLR(:,1),'b-',k,errorLR(:,2),'b--',k,errorGTH(:,1),'r-',k,errorGTH(:,2),'r--');
xlabel('Iteration k');ylabel('Relative error');
legend('LRQBD: E_c', 'LRQBD: E_n','accLRQBD: E_c', 'accLRQBD: E_n');
xlim([0 maxit]);


function [X, Y] = sda(M, n1, opts)
% structured doubling algorithm to solve M-matrix algebraic Riccati equations
%
% [X, Y] = sda(M, I)
%
% compute the solution to XDX - AX - XB + C = 0, where
% M = [B -D; -C A] is a regular M-matrix.
%
% The first block has size n1 x n1
%
% This is SDA/ADDA from Chapter 6 in [Huang, Li, Lin, "Structure-preserving
% doubling algorithms for nonlinear matrix equations"]

arguments
    M
    n1 = floor(length(M)/2)
    opts.safety_factor {mustBeGreaterThanOrEqual(opts.safety_factor, 1)} = 1.2
    opts.maxit = 100
end

if M.left
    error('still to implement by transposing');
end

n = size(M, 1);
d = diag(M);

alphaopt = 1 / max(d(n1+1:end));
betaopt = 1 / max(d(1:n1));

alpha = alphaopt / opts.safety_factor;
beta = betaopt / opts.safety_factor;

M1 = eye(n) + [repmat(alpha, 1, n1) repmat(beta, 1, n-n1)] .* M;
M2 = eye(n) - [repmat(beta, 1, n1) repmat(alpha, 1, n-n1)] .* double(M); % cancellation-safe if opts.safety_factor > 1

% compute M1 \ M2 and M1 \ v with one inversion
EXYFw = M1 \ [M2 (alpha+beta)*M.v];
EXYF = EXYFw(:, 1:end-1); w = EXYFw(:, end);
E = EXYF(1:n1, 1:n1);
X = EXYF(n1+1:end, 1:n1);
Y = EXYF(1:n1, n1+1:end);
F = EXYF(n1+1:end, n1+1:end);

u1 = M.u(1:n1);
u2 = M.u(n1+1:end);
w1 = w(1:n1);
w2 = w(n1+1:end);

k = 0;
while true
    Xold = X;
    Yold = Y;
    % [norm(w1 - (u1 - E*u1 - Y*u2)) ...
    %  norm(w2 -(u2 - X*u1 - F*u2))]  % should be small

    %
    % prepare triplets
    %
    v1 = w1 + E*u1 + Y*(F*u2 + w2);
    v2 = w2 + F*u2 + X*(E*u1 + w1);
    ImYX = mmatrix(-Y*X, u1, v1);
    ImXY = mmatrix(-X*Y, u2, v2);
    % [norm(double(ImYX) - (eye(n1)-Y*X)) ...
    %  norm(double(ImXY) - (eye(n-n1)-X*Y))] % should be small
    
    %
    % main iteration formulas
    %
    EYX = E / ImYX;
    FXY = F / ImXY;
    [E, X, Y, F] = deal(EYX*E, X+FXY*X*E, Y+EYX*Y*F, FXY*F);
    k = k + 1;

    %
    % stopping criterion
    % 

    % since norm(E)*norm(F) tends to 0, the eventually X and Y
    % stop changing in machine precision
    % this componentwise stopping criterion ensures even the small entries
    % have converged
    if all(Xold == X, 'all')
        break
    end
    if nargout > 1 && all(Yold == Y, 'all')
        break
    end
    if k >= opts.maxit
        error('Maximum number of iterations reached')
    end
    
    %
    % prepare w for the next iteration
    %
    [w1, w2] = deal(w1 + EYX*(w1 + Yold*w2), w2 + FXY*(Xold*w1 + w2));
end

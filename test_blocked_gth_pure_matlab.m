clear; clc;

n  = 16;    % Matrix size
nb = 8;     % Block size for blocked algorithms

% Generate test matrices
[Ml, vl, wl] = randtriplet(n, left=true);  % Left case:  v'*M = w'
[Mr, vr, wr] = randtriplet(n);              % Right case: M*v = w

% Run all factorizations:
% left case-----------------------------------------------------------
% gth which corresponds to dmltf2
[L,U] = gthlu(Ml, vl, wl, 'left');
err_gth_left = norm(L*U - Ml, 'fro');

%Recursive blocked gth which corresponds to dmltrf2
[L,U] = gthlu_recursive(Ml, vl, wl, 'left');
err_rgth_left = norm(L*U - Ml, 'fro');

%Blocked gth which corresponds to dmltrf
[L,U] = blocked_gthlu_left(Ml, vl, wl, nb);
err_bgth_left = norm(L*U - Ml, 'fro');

%Recursive blocked gth which corresponds to dmltrff
[L,U] = rec_blocked_gthlu_left(Ml, vl, wl, nb);
err_rbgth_left = norm(L*U - Ml, 'fro');
% right case-----------------------------------------------------------
% gth which corresponds to dmmtf2
[L,U] = gthlu(Mr, vr, wr);
err_gth_right = norm(L*U - Mr, 'fro');

%Recursive blocked gth which corresponds to dmmtrf2
[L,U] = gthlu_recursive(Mr, vr, wr);
err_rgth_right = norm(L*U - Mr, 'fro');

%Blocked gth which corresponds to dmmtrf
[L,U] = blocked_gthlu_right(Mr, vr, wr, nb);
err_bgth_right = norm(L*U - Mr, 'fro');

%Recursive blocked gth which corresponds to dmmtrff
[L,U] = rec_blocked_gthlu_right(Mr, vr, wr, nb);
err_rbgth_right = norm(L*U - Mr, 'fro');

% Display results
fprintf('Errors\n');
fprintf('----------------------------------------------\n');
fprintf('GTH(Left)                     : %e\n',err_gth_left);
fprintf('Recursive  GTH (Left)         : %e\n',err_rgth_left);
fprintf('Blocked GTH (Left)            : %e\n',err_bgth_left);
fprintf('Recursive  blocked GTH (Left) : %e\n',err_rbgth_left);
fprintf('----------------------------------------------\n');
fprintf('GTH(Right)                     : %e\n',err_gth_right);
fprintf('Recursive  GTH (Right)         : %e\n',err_rgth_right);
fprintf('Blocked GTH (Right)            : %e\n',err_bgth_right);
fprintf('Recursive  blocked GTH (Right) : %e\n',err_rbgth_right);
%--------------------------------------------------------------------
% Blocked GTH (left case)
function [L,U] = blocked_gthlu_left(M, v, w, nb)
n = size(M,1);
L = zeros(n);
U = eye(n);
w = w(:)';
for k = 1:nb:n
    kb = min(nb, n-k+1); I = k:k+kb-1; J = k+kb:n;
    if isempty(J)
        w_hat = w(I);
    else
        w_hat = w(I) - v(J) * M(J,I);
    end
    [L11,U11,w_out] = gthlu(M(I,I), v(I), w_hat, 'left');
    L(I,I) = L11;
    U(I,I) = U11;
    if isempty(J), continue; end
    U(I,J) = L11 \ M(I,J);
    L(J,I) = M(J,I) / U11;
    M(J,J) = M(J,J) - L(J,I) * U(I,J);
    % w_scaled=[w(I) / U11 ; w_out  + v(J) * L(J,I)]%
    % w_scaled shows that w_out is not equal to w(I) / U11. We can use w_out+something
    % which contains cancellations or we can solve a triangular system
    w(J) = w(J) - (w(I) / U11) * U(I,J);
end
end
%--------------------------------------------------------------------
% Recursive blocked GTH (left case)
function [L,U] = rec_blocked_gthlu_left(M, v, w, nb)
n = size(M,1);
L = zeros(n);
U = eye(n);
w = w(:)';
for k = 1:nb:n
    kb = min(nb, n-k+1); I = k:k+kb-1; J = k+kb:n;
    if isempty(J)
        w_hat = w(I);
    else
        w_hat = w(I) - v(J) * M(J,I);
    end
    [L11,U11,w_out] = gthlu_recursive(M(I,I), v(I), w_hat, 'left');
    L(I,I) = L11;
    U(I,I) = U11;
    if isempty(J), continue; end
    U(I,J) = L11 \ M(I,J);
    L(J,I) = M(J,I) / U11;
    M(J,J) = M(J,J) - L(J,I) * U(I,J);
    % w_scaled=[w(I) / U11 ; w_out  + v(J) * L(J,I)]%
    % w_scaled shows that w_out is not equal to w(I) / U11. We can use w_out+something
    % which contains cancellations or we can solve a triangular system
    w(J) = w(J) - (w(I) / U11) * U(I,J);
end
end
%--------------------------------------------------------------------
% Blocked GTH (right case)
function [L,U] = blocked_gthlu_right(M, v, w, nb)
n = size(M,1);
L = eye(n);
U = zeros(n);
for k = 1:nb:n
    kb = min(nb, n-k+1);I = k:k+kb-1; J = k+kb:n;
    if isempty(J)
        w_hat = w(I);
    else
        w_hat = w(I) - M(I,J) * v(J);
    end
    [L11,U11,w_out] = gthlu(M(I,I), v(I), w_hat);
    L(I,I) = L11;
    U(I,I) = U11;
    if isempty(J), continue; end
    U(I,J) = L11 \ M(I,J);
    L(J,I) = M(J,I) / U11;
    M(J,J) = M(J,J) - L(J,I) * U(I,J);
    % w_scaled=[L11 \ w(I) , w_out + U(I,J)*v(J)]
    % w_scaled shows that w_out is not equal to L11 \ w(I). We can use w_out+something
    % which contains cancellations or we can solve a triangular system
    w(J) = w(J) - L(J,I) * (L11 \ w(I));
end
end
%--------------------------------------------------------------------
% Recursive blocked GTH  (right case)
function [L,U] = rec_blocked_gthlu_right(M, v, w, nb)
n = size(M,1);
L = eye(n);
U = zeros(n);
for k = 1:nb:n
    kb = min(nb, n-k+1);I = k:k+kb-1; J = k+kb:n;
    [L11,U11,w_out] = gthlu_recursive(M(I,I), v(I),w(I));
    L(I,I) = L11;
    U(I,I) = U11;
    if isempty(J), continue; end
    U(I,J) = L11 \ M(I,J);
    L(J,I) = M(J,I) / U11;
    M(J,J) = M(J,J) - L(J,I) * U(I,J);
    % w_scaled=[L11 \ w(I) , w_out] % these two are equal, so I will think if we can modify dmmtrff based on this
    w(J) = w(J) - L(J,I) * (L11 \ w(I));
end
end
%--------------------------------------------------------------------

function [L, U, w] = gthlu_recursive(M, v, w, isleft)
% Computes an LU factorization of an M-matrix with the GTH trick
%
% [L, U, ww] = gthlu_recursive(M, v, w)
% [L, U, ww] = gthlu_recursive(M, v, w, 'left')
%
% given an M-matrix M and vectors v > 0, w >= 0 such that
%
%   right case:  M*v = w
%   left case:   v'*M = w'
%
% computes a componentwise accurate LU factorization.
%
% The third output ww is:
%   right case: ww = L \ w
%   left case:  ww = w / U

[m, n] = size(M);
assert(m >= n);
assert(length(v)==m);
assert(length(w)==n);

if n == 0
    L = [];
    U = [];
    return
end

if not(exist('isleft','var'))
    isleft = 'right';
end

switch isleft
    case 'left'
        if n==1 % base case
            M(1,1) = (w(1) - v(1,2:end)*M(2:end,1)) / v(1);
            L = [1; M(2:end,1) / M(1,1)];
            U = M(1,1);
            w = w / M(1,1);
            return
        end
        L = eye(m, n);
        U = zeros(n,n);
        n2 = floor(n/2);
        % factor first block
        [L(:,1:n2), U(1:n2,1:n2), w(1,1:n2)] =gthlu_recursive(M(:,1:n2), v, w(1,1:n2), 'left');
        % compute U12
        U(1:n2,n2+1:end) = L(1:n2,1:n2) \ M(1:n2,n2+1:end);
        % update w
        w(1,n2+1:end) = w(1,n2+1:end) - w(1,1:n2) * U(1:n2,n2+1:end);
        % Schur complement
        M(n2+1:end,n2+1:end) = M(n2+1:end,n2+1:end)- L(n2+1:end,1:n2) * U(1:n2,n2+1:end);
        % recurse on second block
        [L(n2+1:end,n2+1:end), U(n2+1:end,n2+1:end), w(1,n2+1:end)] =gthlu_recursive(M(n2+1:end,n2+1:end),v(1,n2+1:end),w(1,n2+1:end),'left');
    case 'right'
        if n == 1 % base case
            U = M(1,1);
            L = [1; M(2:end,1)/U];
            w(2:end) = w(2:end) - L(2:end)*w(1);
            return
        end
        L = eye(m,n);
        U = zeros(n,n);
        n1 = floor(n/2);
        % recursive factorization of A11)
        [L(1:n1,1:n1), U(1:n1,1:n1), w(1:n1)] = gthlu_recursive(M(1:n1,1:n1), v(1:n1), w(1:n1));
        % U12 = L11^{-1} * A12
        U(1:n1,n1+1:n) = L(1:n1,1:n1) \ M(1:n1,n1+1:n);
        % L21 = A21 * U11^{-1}
        L(n1+1:m,1:n1) = M(n1+1:m,1:n1) / U(1:n1,1:n1);
        % A22 = A22 - L21*U12
        M(n1+1:m,n1+1:n) = M(n1+1:m,n1+1:n) - L(n1+1:m,1:n1) * U(1:n1,n1+1:n);
        % w2 = w2 - L21*w1
        w(n1+1:m) = w(n1+1:m) - L(n1+1:m,1:n1) * w(1:n1);
        %  recursive factorization of A22
        [L(n1+1:m,n1+1:n), U(n1+1:n,n1+1:n), w(n1+1:m)] = gthlu_recursive(M(n1+1:m,n1+1:n), v(n1+1:n), w(n1+1:m));
    otherwise
        error('Wrong left/right specification')
end
end


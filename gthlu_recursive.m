function [L, U, w]=gthlu_recursive(M, v, w, isleft)
%computes an LU factorization of an M-matrix with the GTH trick
%
% [L, U, ww] = gthlu(M, v, w)
% [L, U, ww] = gthlu(M, v, w, 'left') %use v'*M=w' instead
%
% given an M-matrix M and vectors v > 0, w >= 0 such that M*v = w, computes 
% a LU factorization in a componentwise accurate way.
% the third optional return parameter is ww = L \ w (or w / U in the 'left'
% case), which is computed almost for free by the algorithm.

[m, n] = size(M);
assert(m >= n);
assert(length(v)==m);
assert(length(w)==n);

%this made some tests fail
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
        %assertVectorsAlmostEqual(v(:)'*M,w(:)','relative',sqrt(eps)*norm(M,1)); %checks that M*v=w holds at least in norm
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
        [L(:,1:n2), U(1:n2,1:n2), w(1,1:n2)] = gthlu_recursive(M(:, 1:n2), v, w(1, 1:n2), 'left');
        U(1:n2, n2+1:end) = L(1:n2,1:n2) \ M(1:n2, n2+1:end);
        w(1, n2+1:end) = w(1, n2+1:end) - w(1,1:n2) * U(1:n2, n2+1:end);
        M(n2+1:end,n2+1:end) = M(n2+1:end,n2+1:end) - L(n2+1:end,1:n2) * U(1:n2, n2+1:end);
        [L(n2+1:end, n2+1:end), U(n2+1:end, n2+1:end), w(1,n2+1:end)] = gthlu_recursive(M(n2+1:end, n2+1:end), v(1, n2+1:end), w(1, n2+1:end), 'left');
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
        [L(1:n1,1:n1), U(1:n1,1:n1), w(1:n1)] = gthlu_recursive(M(1:n1,1:n1), v(1:n1), w(1:n1));
        U(1:n1,n1+1:n) = L(1:n1,1:n1) \ M(1:n1,n1+1:n);
        L(n1+1:m,1:n1) = M(n1+1:m,1:n1) / U(1:n1,1:n1);
        M(n1+1:m,n1+1:n) = M(n1+1:m,n1+1:n) - L(n1+1:m,1:n1) * U(1:n1,n1+1:n);
        w(n1+1:m) = w(n1+1:m) - L(n1+1:m,1:n1) * w(1:n1);
        [L(n1+1:m,n1+1:n), U(n1+1:n,n1+1:n), w(n1+1:m)] = gthlu_recursive(M(n1+1:m,n1+1:n), v(n1+1:n), w(n1+1:m));
    otherwise
        error('Wrong left/right specification');

end        

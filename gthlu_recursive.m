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
            L = M / M(1,1);
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
        error('not implemented yet')
        % L=eye(n);
        % U = M;
        % %assertVectorsAlmostEqual(U*v,w,'relative',sqrt(eps)*norm(M,1)); %checks that M*v=w holds at least in norm
        % for k=1:n-1
        %     pivot = (w(k) - U(k,k+1:end)*v(k+1:end)) / v(k);
        %     U(k,k) = pivot;
        %     L(k+1:end,k) = U(k+1:end,k) / pivot;
        %     U(k+1:end,k) = 0;
        %     w(k+1:end) = w(k+1:end) - L(k+1:end,k)*w(k);
        %     U(k+1:end,k+1:end) = U(k+1:end,k+1:end) - L(k+1:end,k) * U(k,k+1:end); %will write inaccurate values on the diagonal, but we don't care since we'll never use them
        % end
        % U(end,end) = w(end) / v(end);

    otherwise
        error('Wrong left/right specification');

end        

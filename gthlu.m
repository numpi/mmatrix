function [L, U, w]=gthlu(M, v, w, isleft)
%computes an LU factorization of an M-matrix with the GTH trick
%
% [L, U, ww] = gthlu(M, v, w)
% [L, U, ww] = gthlu(M, v, w, 'left') %use v'*M=w' instead
%
% given an M-matrix M and vectors v > 0, w >= 0 such that M*v = w, computes 
% a LU factorization in a componentwise accurate way.
% the third optional return parameter is ww = L \ w (or w / U in the 'left'
% case), which is computed almost for free by the algorithm.

n=length(M);
assert(all(size(M)==[n,n]));
assert(length(v)==n);
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
        L = eye(n); 
        U = M;
        for k = 1:n-1
            pivot = (w(1,k) - v(1,k+1:end) * U(k+1:end, k)) / v(k);
            U(k,k) = pivot;
            L(k+1:end,k) = U(k+1:end,k) / pivot;
            U(k+1:end,k) = 0;
            w(1,k+1:end) = w(1,k+1:end) - w(1,k)/pivot * U(k,k+1:end);
            U(k+1:end,k+1:end) = U(k+1:end,k+1:end) - L(k+1:end,k) * U(k,k+1:end);
        end
        U(end,end) = w(end) / v(end);
        if nargout > 2
            w = w / diag(diag(U));
        end

    case 'right' 
        L=eye(n);
        U = M;
        %assertVectorsAlmostEqual(U*v,w,'relative',sqrt(eps)*norm(M,1)); %checks that M*v=w holds at least in norm
        for k=1:n-1
            pivot = (w(k) - U(k,k+1:end)*v(k+1:end)) / v(k);
            U(k,k) = pivot;
            L(k+1:end,k) = U(k+1:end,k) / pivot;
            U(k+1:end,k) = 0;
            w(k+1:end) = w(k+1:end) - L(k+1:end,k)*w(k);
            U(k+1:end,k+1:end) = U(k+1:end,k+1:end) - L(k+1:end,k) * U(k,k+1:end); %will write inaccurate values on the diagonal, but we don't care since we'll never use them
        end
        U(end,end) = w(end) / v(end);

    otherwise
        error('Wrong left/right specification');

end        

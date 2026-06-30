function [L, U, v]=lu_matlab(M, u, v, opts)
% LU factorization of an M-matrix with the GTH trick (pure-Matlab implementation)
%
% [L, U, w] = lu_matlab(M, u, v)
% [L, U, w] = lu_matlab(M, u, v, left=true) %use u'*M=v' instead
%
% given an M-matrix M and vectors u > 0, v >= 0 such that M*u = v, computes 
% a LU factorization in a componentwise accurate way.
% the third optional return parameter is w = L \ v (or v / U in the 'left'
% case), which is computed almost for free by the algorithm.

arguments
    M (:, :)  % an M-matrix M. The diagonal is recomputed, and assumed to be accurate after construction.
    u (:, :) {mustBePositive} % positive vector	 
    v (:, :) {mustBeNonnegative} % nonnegative vector  
    opts.left (1,1) logical = false
end

n=length(M);
assert(all(size(M)==[n,n]));
assert(length(u)==n);
assert(length(v)==n);

%this made some tests fail
if n == 0
    L = [];
    U = [];
    return
end

if opts.left
    L = eye(n); 
    U = M;
    for k = 1:n-1
        pivot = (v(1,k) - u(1,k+1:end) * U(k+1:end, k)) / u(k);
        U(k,k) = pivot;
        L(k+1:end,k) = U(k+1:end,k) / pivot;
        U(k+1:end,k) = 0;
        v(1,k+1:end) = v(1,k+1:end) - v(1,k)/pivot * U(k,k+1:end);
        U(k+1:end,k+1:end) = U(k+1:end,k+1:end) - L(k+1:end,k) * U(k,k+1:end); %will write inaccurate values on the diagonal, but we don't care since we'll never use them
    end
    U(end,end) = v(end) / u(end);
    if nargout > 2
        v = v / diag(diag(U));
    end
else
    L=eye(n);
    U = M;
    for k=1:n-1
        pivot = (v(k) - U(k,k+1:end)*u(k+1:end)) / u(k);
        U(k,k) = pivot;
        L(k+1:end,k) = U(k+1:end,k) / pivot;
        U(k+1:end,k) = 0;
        v(k+1:end) = v(k+1:end) - L(k+1:end,k)*v(k);
        U(k+1:end,k+1:end) = U(k+1:end,k+1:end) - L(k+1:end,k) * U(k,k+1:end); %will write inaccurate values on the diagonal, but we don't care since we'll never use them
    end
    U(end,end) = v(end) / u(end);
end        

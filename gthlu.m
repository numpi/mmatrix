function [L,U]=gthlu(M,v,w,isleft)
%computes an LU factorization of an M-matrix with the GTH trick
%
% [L,U]=gthlu(M,v,w)
% [L,U]=gthlu(M,v,w,'left') %use v'*M=w' instead
%
% given an M-matrix M and vectors v,w >=0 such that Mv=w, computes L,U
% in a componentwise accurate way.
%
% if 'left', then U has unit diagonal, not L (hack).

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

v = v(:);
w = w(:);

switch isleft
    case 'left'
        %assertVectorsAlmostEqual(v(:)'*M,w(:)','relative',sqrt(eps)*norm(M,1)); %checks that M*v=w holds at least in norm
        L = M;
        U = eye(n);
        v = v'; w = w';
        for k=1:n-1            
            pivot=(w(k)-v(k+1:end)*L(k+1:end,k))/v(k);
            L(k,k)=pivot;
            U(k,k+1:end)=L(k,k+1:end)/pivot;
            L(k,k+1:end)=0;
            L(k+1:end,k+1:end) = L(k+1:end,k+1:end) - L(k+1:end,k)*U(k,k+1:end); %will write inaccurate values on the diagonal, but we don't care since we'll never use them
            w(k+1:end) = w(k+1:end) - w(k)*U(k,k+1:end);
        end
        L(end,end)=w(end)/v(end);
    case 'right' 
        L=eye(n);
        U = M;
        %assertVectorsAlmostEqual(U*v,w,'relative',sqrt(eps)*norm(M,1)); %checks that M*v=w holds at least in norm
        for k=1:n-1
            pivot=(w(k)-U(k,k+1:end)*v(k+1:end))/v(k);
            U(k,k)=pivot;
            L(k+1:end,k)=U(k+1:end,k)/pivot;
            U(k+1:end,k)=0;
            w(k+1:end)=w(k+1:end)-L(k+1:end,k)*w(k);
            U(k+1:end,k+1:end)=U(k+1:end,k+1:end)-L(k+1:end,k)*U(k,k+1:end); %will write inaccurate values on the diagonal, but we don't care since we'll never use them
        end
        U(end,end)=w(end)/v(end);

    otherwise
        error('Wrong left/right specification');
        
        
end        

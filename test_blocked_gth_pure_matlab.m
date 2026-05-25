clear; clc;

n  = 16;    % Matrix size
nb = 8;     % Block size for blocked algorithms

% Generate test matrices
Ml = gettriplet('rand', size=n, left=true);
Mr = gettriplet('rand', size=n, left=false);

% Run all factorizations:
% left case-----------------------------------------------------------
% gth which corresponds to dmltf2
[L,U] = gthlu(Ml.M, Ml.u, Ml.v, 'left');
err_gth_left = norm(L*U - double(Ml), 'fro');

%Recursive blocked gth which corresponds to dmltrf2
[L,U] = gthlu_recursive(Ml.M, Ml.u, Ml.v, 'left');
err_rgth_left = norm(L*U - double(Ml), 'fro');

%Blocked gth which corresponds to dmltrf
[L,U] = blocked_gthlu_left(Ml.M, Ml.u, Ml.v, nb);
err_bgth_left = norm(L*U - double(Ml), 'fro');

%Recursive blocked gth which corresponds to dmltrff
[L,U] = rec_blocked_gthlu_left(Ml.M, Ml.u, Ml.v, nb);
err_rbgth_left = norm(L*U - double(Ml), 'fro');
% right case-----------------------------------------------------------
% gth which corresponds to dmmtf2
[L,U] = gthlu(Mr.M, Mr.u, Mr.v);
err_gth_right = norm(L*U - double(Mr), 'fro');

%Recursive blocked gth which corresponds to dmmtrf2
[L,U] = gthlu_recursive(Mr.M, Mr.u, Mr.v);
err_rgth_right = norm(L*U - double(Mr), 'fro');

%Blocked gth which corresponds to dmmtrf
[L,U] = blocked_gthlu_right(Mr.M, Mr.u, Mr.v, nb);
err_bgth_right = norm(L*U - double(Mr), 'fro');

%Recursive blocked gth which corresponds to dmmtrff
[L,U] = rec_blocked_gthlu_right(Mr.M, Mr.u, Mr.v, nb);
err_rbgth_right = norm(L*U - double(Mr), 'fro');

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
    if isempty(J)
        w_hat = w(I);
    else
        w_hat = w(I) - M(I,J) * v(J);
    end
    [L11, U11, w_out] = gthlu_recursive(M(I,I), v(I),w_hat);
    L(I,I) = L11;
    U(I,I) = U11;
    if isempty(J), continue; end
    U(I,J) = L11 \ M(I,J);
    L(J,I) = M(J,I) / U11;
    M(J,J) = M(J,J) - L(J,I) * U(I,J);
    % w_scaled=[L11 \ w(I) , w_out]
    w(J) = w(J) - L(J,I) * (L11 \ w(I));
end
end
%--------------------------------------------------------------------
 

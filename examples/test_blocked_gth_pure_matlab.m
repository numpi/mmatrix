% This file tests the LU variants (trf, trf2, tf2) in a pure-Matlab implementation

n  = 16;    % Matrix size
nb = 8;     % Block size for blocked algorithms

% Generate test matrices
Ml = gettriplet('rand', size=n, left=true, seed=0);
Mr = gettriplet('rand', size=n, left=false, seed=1);

% Run all factorizations:
% left case-----------------------------------------------------------
% gth which corresponds to dmltf2
[L,U] = lu_matlab(Ml.M, Ml.u, Ml.v, left=true);
err_gth_left = norm(L*U - double(Ml), 'fro');

%Recursive blocked gth which corresponds to dmltrf2
[L,U] = recursive_gthlu_left(Ml.M, Ml.u, Ml.v);
err_rgth_left = norm(L*U - double(Ml), 'fro');

%Blocked gth which corresponds to dmltrf
[L,U] = blocked_gthlu_left(Ml.M, Ml.u, Ml.v, nb);
err_bgth_left = norm(L*U - double(Ml), 'fro');

%Recursive blocked gth which corresponds to dmltrff
[L,U] = rec_blocked_gthlu_left(Ml.M, Ml.u, Ml.v, nb);
err_rbgth_left = norm(L*U - double(Ml), 'fro');
% right case-----------------------------------------------------------
% gth which corresponds to dmmtf2
[L,U] = lu_matlab(Mr.M, Mr.u, Mr.v);
err_gth_right = norm(L*U - double(Mr), 'fro');

%Recursive blocked gth which corresponds to dmmtrf2
[L,U] = recursive_gthlu_right(Mr.M, Mr.u, Mr.v);
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
% Recursive GTH (left case)
function [L, U, w] = recursive_gthlu_left(M, v, w)
    [m, n] = size(M);
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
    [L(:,1:n2), U(1:n2,1:n2), w(1,1:n2)] = recursive_gthlu_left(M(:, 1:n2), v, w(1, 1:n2));
    U(1:n2, n2+1:end) = L(1:n2,1:n2) \ M(1:n2, n2+1:end);
    w(1, n2+1:end) = w(1, n2+1:end) - w(1,1:n2) * U(1:n2, n2+1:end);
    M(n2+1:end,n2+1:end) = M(n2+1:end,n2+1:end) - L(n2+1:end,1:n2) * U(1:n2, n2+1:end);
    [L(n2+1:end, n2+1:end), U(n2+1:end, n2+1:end), w(1,n2+1:end)] = recursive_gthlu_left(M(n2+1:end, n2+1:end), v(1, n2+1:end), w(1, n2+1:end));
end
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
        [L11,U11] = lu_matlab(M(I,I), v(I), w_hat, left=true);
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
        [L11,U11] = recursive_gthlu_left(M(I,I), v(I), w_hat);
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
% Recursive GTH (right case)
function [L, U, w] = recursive_gthlu_right(M, v, w)
    [m, n] = size(M);
    if n == 1 % base case
        U = w(1) / v(1);
        L = [1; M(2:end,1) / U];
        return
    end
    L = eye(m,n);
    U = zeros(n,n);
    n1 = floor(n/2);
    % Save w1 
    w1_old = w(1:n1);
    % Update w1 = w1 - A12 * v2
    w(1:n1) = w(1:n1) - M(1:n1,n1+1:n) * v(n1+1:n);
    % factorization for A11
    [L(1:n1,1:n1), U(1:n1,1:n1), w(1:n1)] = recursive_gthlu_right(M(1:n1,1:n1), v(1:n1), w(1:n1));
    w(1:n1) = w1_old;
    % U12 = L11^{-1} * A12
    U(1:n1,n1+1:n) = L(1:n1,1:n1) \ M(1:n1,n1+1:n);
    % L21 = A21 * U11^{-1}
    L(n1+1:m,1:n1) = M(n1+1:m,1:n1) / U(1:n1,1:n1);
    % Schur complement A22
    M(n1+1:m,n1+1:n) = M(n1+1:m,n1+1:n)- L(n1+1:m,1:n1) * U(1:n1,n1+1:n);
    % w1 = L11^{-1} * w1
    w(1:n1) = L(1:n1,1:n1) \ w(1:n1);
    % w2 = w2 - L21 * w1
    w(n1+1:m) = w(n1+1:m) - L(n1+1:m,1:n1) * w(1:n1);
    [L(n1+1:m,n1+1:n), U(n1+1:n,n1+1:n), w(n1+1:m)] = recursive_gthlu_right(M(n1+1:m,n1+1:n), v(n1+1:n), w(n1+1:m));
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
        [L11,U11] = lu_matlab(M(I,I), v(I), w_hat);
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
        [L11, U11] = recursive_gthlu_right(M(I,I), v(I),w_hat);
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
 

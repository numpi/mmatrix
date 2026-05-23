function [M1, M2] = schur(obj, p, q)
% ------------------------------------------------------------
% [M1, M2] = schur(obj, p, q)
% computes the minor and the Schur complement of M
%   with respect to the indices p  and q in triplet mmatrix form
% M: triplet mmatrix (obj.M)
% p: a vector of indices
% q: the remaining indices (optional)
% M1: the minor as a triplet mmatrix
% M2: the Schur complement as a triplet mmatrix
%
% We assume the matrix M is partitioned as:
%
%   M = [ M11  M12 ]
%       [ M21  M22 ]
% and M u = v (or  u' M = v')
%
% Schur complement:
%  S = M22 - M21 * M11^{-1} * M12
% RHS update:
% for the right triplet: v2 = v2 - M21 * (M11^{-1} v1)
% for the left triplet: v2' = v2' - (v1' * M11^{-1}) * M12
% ------------------------------------------------------------
M = obj.M;
u = obj.u;
v = obj.v;
n = obj.n;

if nargin < 3
    q = setdiff(1:n, p);
end

%  block partition
M11 = M(p,p);
M12 = M(p,q);
M21 = M(q,p);
M22 = M(q,q);
u1 = u(p); u2 = u(q);
v1 = v(p); v2 = v(q);

% Compute the Schur complement
S = M22 - M21 * (M11 \ M12);

if obj.left
    v2 = v2 - (v1 / M11) * M12;
    M1 = mmatrix(M11, u1, v1, left=true);
    M2 = mmatrix(S,   u2, v2, left=true);
else
    v2 = v2 - M21 * (M11 \ v1);
    M1 = mmatrix(M11, u1, v1);
    M2 = mmatrix(S,   u2, v2);
end
end

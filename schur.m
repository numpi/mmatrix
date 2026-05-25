function S = schur(obj, p)
% Schur complement of M with respect to the indices p
%
% S = schur(M, p)
%
% M: mmatrix object
% p: vector of variable indices to remove from M (the "1" block below)
% S: the Schur complement as a triplet mmatrix
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
u = obj.u;
v = obj.v;
n = obj.n;

% create vector of complementary indices
% this syntax works whether p is logical or a list of indices
q = true(1, n);
q(p) = false;

%  block partition
M11 = obj(p,p);
M12 = obj(p,q);
M21 = obj(q,p);
M22 = obj(q,q);
u2 = u(q);
v1 = v(p); v2 = v(q);

% Compute the Schur complement
% we only care about the off-diagonal entries
S = double(M22) + M21 * (M11 \ (-M12));

if obj.left
    v2 = v2 - (v1 / M11) * M12;
    S = mmatrix(S, u2, v2, left=true);
else
    v2 = v2 - M21 * (M11 \ v1);
    S = mmatrix(S, u2, v2);
end
end

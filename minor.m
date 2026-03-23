function [A1, A2] = minor(A, p, q)
% [A1, A2] = minor(A, p, q)
% computes the minor and the Schur complement of A
%   with respect to the indices p  and q in triplet mmatrix form
% A: triplet mmatrix
% p: a vector of indices
% q: the remaining indices (optional)
% A1: the minor as a triplet mmatrix
% A2: the Schur complement as a triplet mmatrix

off = A.off; u = A.u; v = A.v;
m = length(p); n = A.n;

if nargin < 3
    % bugged
    pp = sort(p);
    ct = 1; q = [];
    for ell = 1:n
        ell, pp(ct), pause
        if ne(pp(ct), ell)
            q = [q ell], pause
        else
            ct = ct+1;
        end
    end
end

off1 = off(p, p);
u1 = u(p);
v1 = v(p) - off(p, q)*u(q);
A1 = mmatrix(off1, u1, v1);

off2 = off(q, q) - off(q, p)*(A1\off(p, q));
u2 = u(q);
v2 = v(q) - off(q, p)*(A1\v(p));
A2 = mmatrix(off2, u2, v2);

% U = A_{22} - A_{21}*inv(A_{11})*A_{12}
% A_{11}u_1 + A_{12}u_2 = v_1
% A_{21}u_1 + A_{22}u_2 = v_2
% Uu_2 = A_{22}u_2 - A_{21}*inv(A_{11})*A_{12}u_2
%      = v_2 - A_{21}u_1 - A_{21}*inv(A_{11})*(v_1-A_{11}u_1)
%      = v_2 - A_{21}u_1 - A_{21}*(inv(A_{11})*v_1-u_1)
%      = v_2 - A_{21}u_1 - A_{21}*inv(A_{11})*v_1 + A_{21}u_1
%      = v_2 - A_{21}*inv(A_{11})*v_1 >= 0

end

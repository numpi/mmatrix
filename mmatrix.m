classdef mmatrix
    % mmatrix - a MATLAB class for triplet M-matrices
    %  an n x n M-matrix A is represented through (off, u, v)
    %  off is the matrix of the off-diagonal of A that is
    %  the nonpositive n x n matrix with zero diagonal and whose
    %  off-diagonal coincide with the off-diagonal of A
    %  u is a positive vector, v is a nonegative vector such that A*u = v

    % Properties of the class
    properties (SetAccess = public)
        off (:, :) {mustBeNonpositive} % off-diagonal of A
        u (:, 1) {mustBePositive} % positive vector	such that Au = v
        v (:, 1) {mustBeNonnegative} % nonnegative vector such that Au = v
        left (1,1) logical % true = left, false = right
    end
    properties (Dependent)
        n % size of the triplet mmatrix
        d % diagonal of the triplet mmatrix
    end

    methods % getter, setter, disp

        function n = get.n(A)
            % get the size of the triplet
            n = size(A.off, 1);
        end

        function d = get.d(A)
            % get the diagonal of the triplet mmatrix
            d = (A.v - A.off*A.u)./A.u;
        end

    end

    methods % utilities

        function A = mmatrix(off, u, v, opts)
            % A = mmatrix(off, u, v, ops)
            % Constructor for triplet mmatrix
            % off: off-diagonal elements of the matrix A in a matrix
            % u, v: positive and nonnegative vectors such that A*u = v
            arguments
                off (:,:)
                u (:,1) {mustBePositive}
                v (:,1) {mustBeNonnegative}
                opts.left (1,1) logical = false
            end
            [n, m] = size(off);

            if ne(n, m)
                error('offdiag(A) must be a square matrix');
            end

            if nargin < 4
                opts.left = false;
            end

            A.left = opts.left;


            % if ne(size(u,1), n) || ne(size(v,1), n)
            %     error('u and v should have the same size as offdiag(A)');
            % end

            if A.left
                % left triplet
                if ne(size(u,2), n) || ne(size(v,2), n)
                    error('u and v should have the same size as offdiag(A)');
                end
            else
                % right triplet
                if ne(size(u,1), n) || ne(size(v,1), n)
                    error('u and v should have the same size as offdiag(A)');
                end
            end


            if isequal(diag(off),zeros(n,1))
            else
                warning('A has nonzero diagonal, diagonal elements of A are not used');
                off(1:n+1:n^2) = 0;
            end



            A.off = off;
            A.u = u;
            A.v = v;
        end

        function Ad = double(A)
            % Ad = double(A)
            % makes the triplet mmatrix A a matrix of doubles
            % A: triplet mmatrix
            % Ad: the mmatrix A as a double

            Ad = diag(A.d) + A.off;

        end

        function A2 = vpa(A, f)
            % A2 = vpa(A, f)
            % Makes the triplet mmatrix a vpa matrix
            % A: triplet mmatrix
            % f: digits (optional, default = 40)
            % A2: the triplet mmatrix A as a vpa

            if nargin < 2
                f = 40;
            end

            off = vpa(A.off, f); u = vpa(A.u, f); v = vpa(A.v, f);
            d = (v - off*u)./u;
            A2 = diag(d) + off;

        end

    end

    methods % algorithms

        function [L, U] = lu(A)
            % [L, U] = lu(A)
            % computes the LU factorization of the triplet mmatrix with
            % the GTH algorithm implemented in Fortran
            % A: triplet mmatrix
            % L, U: factors of the LU factorization
            n=A.n;
            if A.left
                U = dmltrff_mex(A.off, A.u, A.v);
            else
                U = dmmtrff_mex(A.off, A.u, A.v);
            end
            L = tril(U,-1) + eye(n);
            U = triu(U);
        end

        function [L, U] = lumatlab(A)
            % [L, U] = lu(A)
            % Forms the LU factorization of A
            % A: triplet mmatrix
            % L, U: factors

            off = A.off; u = A.u; v = A.v;
            n = A.n;

            L = eye(n);
            U = off;

            for ell = 1:n
                U(ell, ell) = (v(ell) - U(ell, ell+1:n)*u(ell+1:n) )/u(ell);
                L(ell+1:n, ell) = U(ell+1:n, ell)/U(ell, ell);
                v(ell+1:n) = v(ell+1:n) - L(ell+1:n, ell)*v(ell);
                U(ell+1:n, ell+1:n) = U(ell+1:n, ell+1:n) - ...
                    L(ell+1:n, ell)*U(ell, ell+1:n);
                U(ell+1:n, ell) = 0;
            end

        end

        function X = mldivide(A, B)
            % X = mldivide(A, B)
            % Computes A\B, accurate if B is nonnegative
            % A: triplet mmatrix
            % B: right hand side
            % X: A\B

            [n1, m] = size(B);

            if ne(n1, A.n)
                error('sizes of A and B are not compatible');
            end

            if sum(sum(B<0)) > 0
                warning('B is no nonnegative, result may not be componentwise accurate');
            end
            if A.left
                X = dmltrff_mexsolve(A.off, A.u, A.v, B);
            else
                X = dmmtrff_mexsolve(A.off, A.u, A.v, B);
            end
            % [L, U] = lu(A);
            %
            % % forward substitution
            % y = zeros(A.n, m);
            % for ell = 1:A.n
            %     y(ell, :) = B(ell, :) - L(ell, 1:ell-1)*y(1:ell-1, :);
            % end
            %
            % X = zeros(A.n, m);
            % % back substitution
            % for ell = A.n:-1:1
            %     X(ell, :) = (y(ell, :) - U(ell, ell+1:A.n)*X(ell+1:A.n, :)) / U(ell,ell);
            % end
        end

        function X = inv(A)
            % X = inv(A)
            % computes the inverse of the triplet mmatrix A
            % A: triplet mmatrix
            % X: inv(A)
            X = A\eye(A.n);
            % X = dmmtri_mex(A.off, A.u, A.v);
        end

        %		function [A1, A2] = minor(A, p, q)

        function v = perron(A)
            % Not defined

        end

        function B = sqrtm(A)
            % Not defined
        end

    end

    methods %

        function C = plus(A, B)
            warning('used double for the sum, results may not be componentwise accurate')
            C = double(A) + double(B);
        end
        function C = minus(A, B)
            warning('used double for the subtraction, results may not be componentwise accurate')
            C = double(A) - double(B);
        end
        function C = mprod(A, B)
            warning('used double for the product, results may not be componentwise accurate')
            C = double(A) * double(B);
        end
        function x = norm(A)
            warning('norm undefined for triplet mmatrix, used double(A)')
            x = norm(double(A));
        end
    end
end

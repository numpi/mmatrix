classdef mmatrix
    % mmatrix - a MATLAB class for triplet M-matrices
    %  an n x n M-matrix M is represented through (M, u, v)
    %  u is a positive vector, v is a nonegative vector such that M*u = v

    % Properties of the class
    properties (SetAccess = public)
        M (:, :)   % an M-matrix M
        u (:, :) {mustBePositive} % positive vector	such that Mu = v
        v (:, :) {mustBeNonnegative} % nonnegative vector such that Mu = v
        left (1,1) logical % true = left, false = right
    end
    properties (Dependent)
        n % size of the triplet mmatrix
    end

    methods % getter, setter, disp

        function n = get.n(obj)
            % get the size of the triplet
            n = size(obj.M, 1);
        end

    end

    methods % utilities

        function obj = mmatrix(M, u, v, opts)
            % obj = mmatrix(M, u, v, ops)
            % Constructor for triplet mmatrix
            % M: an M-matrix
            % u, v: positive and nonnegative vectors such that M*u = v
            arguments
                M (:,:)
                u (:,1) {mustBePositive}
                v (:,1) {mustBeNonnegative}
                opts.left (1,1) logical = false
            end
            [n, m] = size(M);

            if ne(n, m)
                error('M must be a square matrix');
            end

            if nargin < 4
                opts.left = false;
            end
            obj.left = opts.left;
            if obj.left
                % left triplet
                if ne(size(u,2), n) || ne(size(v,2), n)
                    error('u and v should have the same size as M');
                end
            else
                % right triplet
                if ne(size(u,1), n) || ne(size(v,1), n)
                    error('u and v should have the same size as M');
                end
            end
                
            M(1:n+1:n^2) = 0;
            M(1:n+1:n^2) = (v - M*u)./u;
            
            obj.M = M;
            obj.u = u;
            obj.v = v;
        end

        function Ad = double(obj)
            % Ad = double(A)
            % makes the triplet mmatrix A a matrix of doubles
            % A: triplet mmatrix
            % Ad: the mmatrix A as a double
            Ad = obj.M;
        end

        function A2 = vpa(obj, f)
            % A2 = vpa(A, f)
            % Makes the triplet mmatrix a vpa matrix
            % A: triplet mmatrix
            % f: digits (optional, default = 40)
            % A2: the triplet mmatrix A as a vpa
            if nargin < 2
                f = 40;
            end
            M = vpa(obj.M, f); u = vpa(obj.u, f); v = vpa(obj.v, f);
            M(1:n+1:n^2) = 0;
            d = (v - M*u)./u;
            A2 = diag(d) + M;
        end
    end

    methods % algorithms

        function [L, U] = lu(obj)
            % [L, U] = lu(A)
            % computes the LU factorization of the triplet mmatrix with
            % the GTH algorithm implemented in Fortran
            % A: triplet mmatrix
            % L, U: factors of the LU factorization
            if obj.left
                U = dmltrff_mex(obj.M, obj.u, obj.v);
            else
                U = dmmtrff_mex(obj.M, obj.u, obj.v);
            end
            L = tril(U,-1) + eye(obj.n);
            U = triu(U);
        end

        function [L, U] = lumatlab(obj)
            % [L, U] = lu(A)
            % Forms the LU factorization of A
            % A: triplet mmatrix
            % L, U: factors
            M = obj.M; u = obj.u; v = obj.v;
            n = obj.n;
            L = eye(n);
            U = M;
            for ell = 1:n
                U(ell, ell) = (v(ell) - U(ell, ell+1:n)*u(ell+1:n) )/u(ell);
                L(ell+1:n, ell) = U(ell+1:n, ell)/U(ell, ell);
                v(ell+1:n) = v(ell+1:n) - L(ell+1:n, ell)*v(ell);
                U(ell+1:n, ell+1:n) = U(ell+1:n, ell+1:n) - ...
                    L(ell+1:n, ell)*U(ell, ell+1:n);
                U(ell+1:n, ell) = 0;
            end
        end

        function X = mldivide(obj, B)
            % X = mldivide(A, B)
            % Computes A\B, accurate if B is nonnegative
            % A: triplet mmatrix
            % B: right hand side
            % X: A\B
            [n1, m] = size(B);
            if ne(n1, obj.n)
                error('sizes of A and B are not compatible');
            end

            if sum(sum(B<0)) > 0
                warning('B is no nonnegative, result may not be componentwise accurate');
            end
            if obj.left
                X = dmltrff_mexsolve(obj.M, obj.u, obj.v, B);
            else
                X = dmmtrff_mexsolve(obj.M, obj.u, obj.v, B);
            end
        end


        function X = inv(M)
            % X = inv(M)
            % computes the inverse of the triplet mmatrix M
            % M: triplet mmatrix
            % X: inv(M)
            X = M\eye(obj.n);
            % X = dmmtri_mex(obj.M, obj.u, obj.v);
        end
        %		function [A1, A2] = minor(A, p, q)

        function X = mrdivide(obj, B)  
              % [n1, m] = size(B);
  
            if obj.left
                X = dmltrff_solve(obj.M, obj.u, obj.v, B');
            else
                X = dmmtrff_solve(obj.M, obj.u, obj.v, B');
            end
        end
 
        function A = subrefs(M,S)
            

            A = M.M(S)
        end

        function v = perron(M)
            % Not defined
        end

        function B = sqrtm(M)
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

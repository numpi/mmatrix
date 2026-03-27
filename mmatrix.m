classdef mmatrix
    % mmatrix - a MATLAB class for triplet M-matrices
    %  an n x n M-matrix M is represented through (M, u, v)
    %  u is a positive vector, v is a nonegative vector

    % Properties of the class
    properties (SetAccess = public)
        M (:, :)   % an M-matrix M
        u (:, :) {mustBePositive} % positive vector	 
        v (:, :) {mustBeNonnegative} % nonnegative vector  
        left (1,1) logical
        % left = true  --> u' M = v'
        % left = false --> M u = v
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
        %------------------------------------------------------------------

    methods % utilities

        function obj = mmatrix(M, u, v, opts)
            % obj = mmatrix(M, u, v, ops)
            % Constructor for triplet mmatrix
            % M: an M-matrix
            % u, v: positive and nonnegative vectors such that M*u = v(or u' M = v')
            arguments
                M (:,:)
                u (:,:) {mustBePositive}
                v (:,:) {mustBeNonnegative}
                opts.left (1,1) logical = false
            end
            [n, m] = size(M);
            if ne(n, m)
                error('M must be a square matrix');
            end
            obj.left = opts.left;

            if obj.left
                u = u(:)'; v = v(:)';
            else
                u = u(:);  v = v(:);
            end
            if length(u) ~= n || length(v) ~= n
                error('u and v should have length n');
            end
            M(1:n+1:n^2) = 0;
            if obj.left
                M(1:n+1:n^2) = (v - u*M )./u;
            else
                M(1:n+1:n^2) = (v - M*u) ./ u;
            end
            obj.M = M;
            obj.u = u;
            obj.v = v;
        end
        %------------------------------------------------------------------
        function Md = double(obj)
            % Md = double(M)
            % makes the triplet mmatrix M a matrix of doubles
            % M: triplet mmatrix
            % Md: the mmatrix M as a double
            Md = obj.M;
        end
        %------------------------------------------------------------------
        function M2 = vpa(obj, f)
            % M2 = vpa(M, f)
            % Makes the triplet mmatrix a vpa matrix
            % M: triplet mmatrix
            % f: digits (optional, default = 40)
            % M2: the triplet mmatrix M as a vpa
            if nargin < 2
                f = 40;
            end
            n=obj.n;
            M = vpa(obj.M, f);
            u = vpa(obj.u, f);
            v = vpa(obj.v, f);
            M(1:n+1:n^2) = 0;
            if obj.left
                d = (v - u*M )./u;
            else
                d = (v - M*u) ./ u;
            end
            M2 = diag(d) + M;
        end
    end

    methods % algorithms
        %------------------------------------------------------------------
        function [L, U] = lu(obj)
            % [L, U] = lu(M)
            % computes the LU factorization of the triplet mmatrix with
            % the GTH algorithm implemented in Fortran
            % M: triplet mmatrix
            % L, U: factors of the LU factorization
            if obj.left
                U = dmltrff(obj.M, obj.u, obj.v);
            else
                U = dmmtrff(obj.M, obj.u, obj.v);
            end
            L = tril(U,-1) + eye(obj.n);
            U = triu(U);
        end
        %------------------------------------------------------------------
        function [L, U] = lumatlab(obj)
            % [L, U] = lu(M)
            % Forms the LU factorization of M
            % M: triplet mmatrix
            % L, U: factors
            if obj.left
                [L, U] = gthlu(obj.M, obj.u, obj.v,'left');
            else
                [L, U] = gthlu(obj.M, obj.u, obj.v);
            end
        end
        %------------------------------------------------------------------
        function X = mldivide(obj, B)
            % X = mldivide(M, B)
            % Computes M\B, accurate if B is nonnegative
            % M: triplet mmatrix
            % B: right hand side
            % X: M\B
            n1= size(B,1);
            if ne(n1, obj.n)
                error('sizes of M and B are not compatible');
            end
            if sum(sum(B<0)) > 0
                warning('B is no nonnegative, result may not be componentwise accurate');
            end
            if obj.left
                X = dmltrff_solve(obj.M, obj.u, obj.v, B);
            else
                X = dmmtrff_solve(obj.M, obj.u, obj.v, B);
            end
        end
        %------------------------------------------------------------------
        function X = inv(obj)
            % X = inv(M)
            % computes the inverse of the triplet mmatrix M
            % M: triplet mmatrix
            % X = obj.M\eye(obj.n);
            if obj.left
                % Mf = dmltrff(obj.M, obj.u, obj.v);
                % X = dmmtri_mex(Mf, obj.u, obj.v); this works too.
                X = dmltrff_solve(obj.M, obj.u, obj.v, eye(obj.n));
            else
                Mf = dmmtrff(obj.M, obj.u, obj.v); 
                X = dmmtri_mex(Mf, obj.u, obj.v);
                % equivalently, we can solve a linear system with I in the right side,
                % but dmmtri seems to be faster. 
                % X = dmltrff_solve(obj.M, obj.u, obj.v, eye(obj.n)); 
            end
        end
        %------------------------------------------------------------------
        function X = mrdivide(obj, B)
            % X = mldivide(M, B)
            % Computes M\B, accurate if B is nonnegative
            % M: triplet mmatrix
            % B: left hand side
            % X: B/M
            m1 = size(B,2);
            if ne(m1, obj.n)
                error('sizes of M and B are not compatible');
            end
            if sum(sum(B<0)) > 0
                warning('B is no nonnegative, result may not be componentwise accurate');
            end
            if obj.left
                X = dmltrff_lsolve(obj.M, obj.u, obj.v, B);
            else
                X = dmmtrff_lsolve(obj.M, obj.u, obj.v, B);
            end
        end
        %------------------------------------------------------------------
        function A = subrefs(M,S)
            A = M.M(S)
        end
        %------------------------------------------------------------------
        function v = perron(M)
            % Not defined
        end
        %------------------------------------------------------------------
        function B = sqrtm(M)
            % Not defined
        end
    end

    %----------------------------------------------------------------------
    methods %
        function disp(obj)
            fprintf('\nM:\n');
            if obj.n <= 10
                disp(obj.M)
            else
                disp(obj.M(1:10,1:10))
                fprintf('... (truncated)\n');
            end
            fprintf('\nu:\n');
            disp(obj.u.');
            fprintf('v:\n');
            disp(obj.v.');
        end
        %------------------------------------------------------------------
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

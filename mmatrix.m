classdef mmatrix <  matlab.mixin.indexing.RedefinesParen
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
            % Constructor for triplet mmatrix
            %
            % obj = mmatrix(M, u, v)
            % obj = mmatrix(M, u, v, left=true)
            %
            % M: an M-matrix
            % u, v: positive and nonnegative vectors such that M*u = v 
            % (or u' M = v' if left=true)
            %
            % The provided diagonal entries of M are ignored and recomputed from u
            % and v instead.
            %
            % Optionally, v may be computed as M*u or u*M, if not provided,
            % however, this is only componentwise accurate if M is diagonal
            % and will give a warning otherwise:
            %
            % obj = mmatrix(D, u)
            % obj = mmatrix(D, u, left=true)

            arguments
                M (:,:)
                u (:,:) {mustBePositive}
                v (:,:) {mustBeNonnegative} = []
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

            if isempty(v)
                if not(isdiag(M))
                    warning('mmatrix:inaccurate', 'computed v from M, results may not be componentwise accurate')
                end
                if obj.left
                    v = u*M;
                else
                    v = M*u;
                end
                if not(all(v>=0))
                    error('invalid triplet: the product between M and u must be nonnegative')
                end
            end

            if length(u) ~= n || length(v) ~= n
                error('u and v should have length n');
            end
            M(1:n+1:n^2) = 0;
            if not(all(M(:)<=0))
                error('invalid triplet: M cannot have positive off-diagonal entries');
            end
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
            n = obj.n;
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
                U = mmex.dmltrff(obj.M, obj.u, obj.v);
            else
                U = mmex.dmmtrff(obj.M, obj.u, obj.v);
            end
            L = tril(U,-1) + eye(obj.n);
            U = triu(U);
        end
        %------------------------------------------------------------------
               %------------------------------------------------------------------
        function [L, U] = lumatlab(obj)
            % [L, U] = lumatlab(M)
            % Forms the LU factorization of M using pure MATLAB
            % M: triplet mmatrix
            % L, U: factors
                        
            n = obj.n;
            if obj.left
                L = eye(n); 
                U = obj.M;
                u = obj.u;
                v = obj.v;

                for j = 1:n-1
                    pivot = (v(1, j) - u(1, j+1:n) * U(j+1:n, j)) / u(1, j);
                    U(j, j) = pivot;
                    L(j+1:n, j) = U(j+1:n, j) / pivot;
                    U(j+1:n, j) = 0;
                    v(1, j+1:n) = v(1, j+1:n) - v(1, j)/pivot * U(j, j+1:n);
                    U(j+1:n, j+1:n) = U(j+1:n, j+1:n) - ...
                                      L(j+1:n, j) * U(j, j+1:n);
                end
                U(n, n) = v(1, n) / u(1, n);
            else
                L = eye(n);
                U = obj.M;
                u = obj.u;
                v = obj.v;
                for j = 1:n-1
                    pivot = (v(j) - U(j, j+1:n) * u(j+1:n)) / u(j);
                    U(j, j) = pivot;
                    L(j+1:n, j) = U(j+1:n, j) / pivot;
                    U(j+1:n, j) = 0;
                    v(j+1:n) = v(j+1:n) - L(j+1:n, j) * v(j);
                    U(j+1:n, j+1:n) = U(j+1:n, j+1:n) - ...
                                      L(j+1:n, j) * U(j, j+1:n); 
                    % will write inaccurate values on the diagonal, 
                    % but we don't care since we'll never use them
                end
                U(n, n) = v(n) / u(n);
            end
        end
        %------------------------------------------------------------------
        function [L, U] = lumatlabold(obj)
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
            if any(B(:) < 0)
                warning('mmatrix:inaccurate', 'B is not nonnegative, result may not be componentwise accurate');
            end
            if obj.left
                X = mmex.dmltrff_solve(obj.M, obj.u, obj.v, B);
            else
                X = mmex.dmmtrff_solve(obj.M, obj.u, obj.v, B);
            end
        end
        %------------------------------------------------------------------
        function X = inv(obj)
            % X = inv(M)
            % computes the inverse of the triplet mmatrix M
            % M: triplet mmatrix
            % X = obj.M\eye(obj.n);
            % here we can use dmmtri for the both left and right triplet
            % because it does not depent on u and v
            if obj.left
                Mf = mmex.dmltrff(obj.M, obj.u, obj.v);
            else
                Mf = mmex.dmmtrff(obj.M, obj.u, obj.v);
            end
            X = dmmtri_mex(Mf);
        end
        %------------------------------------------------------------------
        function X = mrdivide(B, obj)
            % X = mrdivide(B, M)
            % Computes B/M, accurate if B is nonnegative
            % M: triplet mmatrix
            % B: left hand side
            % X: B/M
            m1 = size(B,2);
            if ne(m1, obj.n)
                error('sizes of M and B are not compatible');
            end
            if any(B(:) < 0)
                warning('mmatrix:inaccurate', 'B is not nonnegative, result may not be componentwise accurate');
            end
            if obj.left
                X = mmex.dmltrff_lsolve(obj.M, obj.u, obj.v, B);
            else
                X = mmex.dmmtrff_lsolve(obj.M, obj.u, obj.v, B);
            end
        end
        %------------------------------------------------------------------
        function X = det(obj)
        % Determinant of a triplet M-matrix
        % X = det(obj) returns det(M) computed from the GTH algorithm
        % Since L is unit lower triangular, det(M) = prod(diag(U)).
            [~, U] = lu(obj);
            X = prod(diag(U));
        end
        %------------------------------------------------------------------
        function [U, S, V] = svd(obj)
            % [U, S, V] = svd(M)
            % Computes SVD of an M-matrix

            % This code implements Algorithm 1 from:
            % J. Demmel and P. Koev,
            % "Accurate SVDs of weakly diagonally dominant M-matrices,"
            % Numerische Mathematik, 98:99–104, 2004.
            % to obtain LDU factorization and Algorithm 3.1 from:
            % J. Demmel, et al. "Computing the singular value decomposition
            % with high relative accuracy," Linear Algebra and its
            % Applications, 299 (1999), 21–80.

            A_off = obj.M;
            A_off(1:obj.n+1:end) = 0;
            s = zeros(obj.n,1);
            for i = 1:obj.n
                s(i) = obj.M(i,i);
                for j = 1:obj.n
                    if j ~= i
                        s(i) = s(i) + obj.M(i,j);
                    end
                end
            end
            [U, sigma, V] = dmmsvd_mex(A_off, s);
            S = diag(sigma);
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function [lambda, u, iter, info] = eigmin(obj, tol, maxit)
        % Smallest eigenvalue of a triplet M-matrix

        % Inputs:
        % obj   : mmatrix object (triplet M-matrix)
        % tol   : convergence tolerance (default 100*eps)
        % maxit : maximum number of iterations (default 100)

        % Outputs:
        % lambda : approximation of smallest eigenvalue
        % u      : corresponding eigenvector  
        % iter   : number of iterations performed
        % info   : structure with error analysis
        
        % This code implements Algorithm 2 from:
        % A. S. Alfa, J. Xue, and Q. Ye, "Accurate Computation of the
        % Smallest Eigenvalue of a Diagonally Dominant M-Matrix,"
        % Mathematics of Computation, Vol. 71, No. 237, pp. 217–236, 2001.

            if nargin < 2 || isempty(tol)
                tol = 100*eps;
            end
            if nargin < 3
                maxit = 100;
            end
            u_s = obj.u;
            lambda = min(obj.v ./ obj.u);
            v_s = obj.v - lambda * obj.u;

            for iter = 1:maxit

                if obj.left
                    w = mmex.dmltrff_lsolve(obj.M, u_s, v_s, u_s);
                else
                    w = mmex.dmmtrff_solve(obj.M, u_s, v_s, u_s);
                end

                ratios = u_s ./ w;
                theta = min(ratios);

                lambda_new = lambda + theta;

                u_new = w / norm(w, inf);
                v_new = u_new .* (ratios - theta);

                if (max(ratios) - min(ratios)) / lambda_new <= tol
                    lambda = lambda_new;
                    u = u_new;
                    break;
                end

                lambda = lambda_new;
                u_s = u_new;
                v_s = v_new;
            end

            if iter == maxit
                warning('mmatrix:eigmin', 'Max iterations reached');
            end
            
            % Error analysis
            eps_machine = eps;
            phi = (2*(obj.n+2)*(obj.n+3)*(2*obj.n+5))/3;
            varphi = 12*obj.n*(phi + 3); 
            psi = 12*(iter+1)*obj.n*varphi + 2*phi + 36*(iter+1)*obj.n + 4; 

            info.tol = tol;
            info.phi = phi;
            info.varphi = varphi;
            info.psi = psi;

            % error bound
            info.bound = tol + psi * eps_machine;
        end
        %------------------------------------------------------------------
        function [B, out] = sqrtm(M, opts)
            % Computes M^{1/2}, the square root of the M-matrix M
            % B = sqrtm(M) 
            % [B, out] = sqrtm(M)
            % B = sqrtm(M, g = 4, maxit = 20, tol = 1d-20)
            % 
            % M: triplet M-matrix
            % opts.gamma: cancellation reducing parameters (default = 4)
            % opts.maxit: maximum number of iterations (default = 20)
            % opts.tol: tolerance in the stopping criterion (default = 1d-20)
            % B: the square root
            % out: output structure, contains
            %      out.it: number of iterations
            %
            % The algorithm is Newton's method written in Cyclic Reduction
            % form and with an iteration on the triplet, from
            %
            % D. A. Bini, B. Iannazzo, B. Meini, J. Meng, "Componentwise
            % accurate computation of the square root of an M-matrix",
            % manuscript, 2026.
            %
            % Warning: the algorithm makes a subtraction in the first step,
            % the parameter gamma reduces the cancellation to subtraction 
            % of the type (gamma - a) with a<=1, an high value of gamma
            % slows down the convergence.

            arguments
                M
                opts.gamma = 4
                opts.maxit = 20
                opts.tol = 1d-20
            end
            gamma = opts.gamma;

            n = M.n;
            gamma = gamma * max(diag(M));
            M = 1/gamma * M;
            W = M.M - eye(n); % not componentwise accurate
            Z = (1/2) * (eye(n) + M);
            p = M.v;
            for ell = 1:opts.maxit
                [L, U] = lumatlab(Z);
                p = p - 1/2 * W * (U\(L\p)); % p yields the triplet vector v of Z
                W = - 1/4 * W * (U\(L\W)); % W is the update on the off-diagonal
                Z = mmatrix(Z.M + 1/2 * W, M.u, p - 1/2 * W * M.u);
                if norm(W) < opts.tol
                   break;
                end
            end
            out.it = ell;
            B = sqrt(gamma) * Z;
        end

    end

    %----------------------------------------------------------------------
    methods %
        function disp(obj)
            fprintf('M-matrix with entries\n');
            if obj.n <= 10
                disp(obj.M)
            else
                disp(obj.M(1:10,1:10))
                fprintf('... (truncated)\n');
            end
            if obj.left
                relation = 'u*M = v';
            else
                relation = 'M*u = v';
            end
            fprintf('satisfying %s with u =\n', relation);
            disp(obj.u);
            fprintf('and v=\n');
            disp(obj.v);
        end
        %------------------------------------------------------------------
        function C = plus(A, B)
            % Matrix sum with accurate paths for M+M, M+Diagonal, Diagonal+M
            if isa(A, 'mmatrix') && isa(B, 'mmatrix') && A.left == B.left && isequal(A.u, B.u)
                C = mmatrix(A.M+B.M, A.u, A.v+B.v, left=A.left);
            elseif isa(A, 'mmatrix') && isa(B, 'double') && isdiag(B) && all(diag(B)>=0)
                C = A + mmatrix(B, A.u, left=A.left);
            elseif isa(A, 'double') && isdiag(A) && all(diag(A)>=0) && isa(B, 'mmatrix')
                C = mmatrix(A, B.u, left=B.left) + B;
            else
                warning('mmatrix:inaccurate', 'converted to double for the sum, results may not be componentwise accurate')
                C = double(A) + double(B);
            end
        end
        function C = minus(A, B)
            warning('mmatrix:inaccurate', 'converted to double for the subtraction, results may not be componentwise accurate')
            C = double(A) - double(B);
        end
        function C = times(A, B)
            if size(A, 1) == 1 && size(A, 2) == 1 && A(1, 1) > 0
                C = mmatrix(A*B.M, B.u, A*B.v);
            elseif size(B, 1) == 1 && size(B, 2) == 1 && B(1, 1) > 0
                C = mmatrix(A.M*B, A.u, A.v*B);
            else
                warning('mmatrix:inaccurate', 'converted to double for the product, results may not be componentwise accurate')
                C = double(A) .* double(B);
            end
        end
        function C = mtimes(A, B)
            if size(A, 1) == 1 && size(A, 2) == 1 && A(1, 1) > 0
                C = mmatrix(A*B.M, B.u, A*B.v);
            elseif size(B, 1) == 1 && size(B, 2) == 1 && B(1, 1) > 0
                C = mmatrix(A.M*B, A.u, A.v*B);
            else
                warning('mmatrix:inaccurate', 'converted to double for the product, results may not be componentwise accurate')
                C = double(A) * double(B);
            end
        end
        function x = norm(A, varargin)
            x = norm(double(A), varargin{:});
        end

    end

    % the following methods are required from matlab.mixin.indexing.RedefinesParen
    methods
        function C = cat(dim, varargin)
            error('Concatenation is unsupported: construct a new mmatrix')
        end
        function e = empty(varargin)
            e = mmatrix([], [], []);
        end
        function sz = size(obj, varargin)
            sz = size(obj.M, varargin{:});
        end
    end

    % the following methods are required from matlab.mixin.indexing.RedefinesParen
    methods (Access=protected)
        function out = parenAssign(varargin)
            error('Modification is unsupported: construct a new mmatrix');
        end
        function parenDelete(varargin)
            error('Modification is unsupported; construct a new mmatrix');
        end        
        function n = parenListLength(obj,indexOp,ctx)
            n = parenListLength(obj.M, indexOp, ctx);
        end
        function result = parenReference(obj, indexOp)
            if length(indexOp.Indices)==2 && isequal(indexOp.Indices(1), indexOp.Indices(2))
                I = indexOp.Indices{1};
                notI = true(size(obj,2), 1);
                notI(I) = false;
                if obj.left
                    result = mmatrix(obj.M(I,I), obj.u(I), obj.v(I) - obj.u(notI)*obj.M(notI,I), left=true);
                else
                    result = mmatrix(obj.M(I,I), obj.u(I), obj.v(I) - obj.M(I,notI)*obj.u(notI), left=false);
                end
            else
                result = obj.M.(indexOp);
            end
        end
    end
end

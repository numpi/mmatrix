classdef mmatrix <  matlab.mixin.indexing.RedefinesParen
    % mmatrix - a MATLAB class for triplet M-matrices
    %  an n x n M-matrix M is represented through (M, u, v)
    %  u is a positive vector, v is a nonegative vector

    % Properties of the class
    properties (SetAccess = public)
        M (:, :)  % an M-matrix M. The diagonal is recomputed, and assumed to be accurate after construction.
        u (:, :) {mustBePositive} % positive vector	 
        v (:, :) {mustBeNonnegative} % nonnegative vector  
        left (1,1) logical
        % left = true  --> u' M = v'
        % left = false --> M u = v
    end
    properties (Dependent)
        n % size of the triplet mmatrix
    end

    methods % getter, setter

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
            % Convert the triplet mmatrix M intoa matrix of doubles
            % M: triplet mmatrix
            % Md: the mmatrix M as a double
            Md = obj.M;
        end
        %------------------------------------------------------------------
        function M2 = vpa(obj, f)
            % M2 = vpa(M, f)
            % Convert the triplet mmatrix into a vpa matrix
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
        function [L, U] = lu(obj, opts)
            % [L, U] = lu(M)
            % compute the LU factorization of the triplet mmatrix with
            % the GTH algorithm implemented in Fortran
            % M: triplet mmatrix
            % L, U: factors of the LU factorization
            arguments
                obj
                opts.use_mex (1,1) logical = true
            end
            if opts.use_mex
                if obj.left
                    U = mmex.dmltrff(obj.M, obj.u, obj.v);
                else
                    U = mmex.dmmtrff(obj.M, obj.u, obj.v);
                end
                L = tril(U,-1) + eye(obj.n);
                U = triu(U);
            else
                [L, U] = lu_matlab(obj.M, obj.u, obj.v, obj.left);
            end
        end
        %------------------------------------------------------------------
        function X = mldivide(obj, B)
            % X = mldivide(M, B)
            % Compute M\B, accurate if B is nonnegative
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
            % Compute the inverse of the triplet mmatrix M
            % M: triplet mmatrix
            % X = obj.M\eye(obj.n);
            % here we can use dmmtri for the both left and right triplet
            % because it does not depent on u and v
            if obj.left
                Mf = mmex.dmltrff(obj.M, obj.u, obj.v);
            else
                Mf = mmex.dmmtrff(obj.M, obj.u, obj.v);
            end
            X = mmex.dmmtri(Mf);
        end
        %------------------------------------------------------------------
        function X = mrdivide(B, obj)
            % X = mrdivide(B, M)
            % Compute B/M, accurately if B is nonnegative
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

            if obj.left
                [V, sigma, U] = mmex.dmmsvd(obj.M', obj.u', obj.v');
            else
                [U, sigma, V] = mmex.dmmsvd(obj.M, obj.u, obj.v);
            end
            
            if nargout == 1
                U = sigma;
            else
                S = diag(sigma);
            end
        end
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

            if all(obj.v == 0)
                lambda = 0;
                u = obj.u;
                iter = 0;
                return
            end

            if nargin < 2 || isempty(tol)
                tol = 100*eps;
            end
            if nargin < 3
                maxit = 100;
            end
            u = obj.u; v = obj.v;
            lambda = min(v./u);
            v = v - lambda * u;
            for iter = 1:maxit
                if obj.left
                    z = mmex.dmltrff_lsolve(obj.M, u, v, u);
                else
                    z = mmex.dmmtrff_solve(obj.M, u, v, u);
                end
                p = u ./ z;
                lambda = lambda + min(p);
                u = z / norm(z, inf);
                v = u .* (p - min(p));
                if (max(p) - min(p)) / lambda <= tol
                    break;
                end
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
            % Compute M^{1/2}, the square root of the M-matrix M
            % B = sqrtm(M) 
            % [B, out] = sqrtm(M)
            % B = sqrtm(M, gamma = 1.5, maxit = 100, tol = 1d-20)
            % 
            % M: triplet M-matrix
            % opts.gamma: cancellation reducing parameters (default = 1.5)
            % opts.maxit: maximum number of iterations (default = 100)
            % opts.tol: tolerance in the stopping criterion (default = 1d-20)
            % B: the square root
            % out: output structure, contains
            %      out.it: number of iterations
            %
            % The algorithm is Newton's method written in Cyclic Reduction
            % form and with an iteration on the triplet, from
            %
            % D. A. Bini, B. Iannazzo, B. Meini, J. Meng, "Component-wise
            % accurate computation of the square root of an M-matrix",
            % arXiv, 2026 (https://arxiv.org/abs/2605.21679).
            %
            % Warning: the algorithm makes a subtraction in the first step,
            % the parameter gamma reduces the cancellation to subtraction 
            % of the type (1 - 1/gamma) with a<=1, an high value of gamma
            % slows down the convergence.

            arguments
                M
                opts.gamma = 1.5
                opts.maxit = 100
                opts.tol = 1d-20
                opts.use_mex (1,1) logical = true
            end

            if M.left
                % calls sqrtm(M') with the same arguments
                nv = namedargs2cell(opts);
                [B, out] = sqrtm(M', nv{:});
                B = B';
                return
            end

            gamma = opts.gamma;

            n = M.n;
            gamma = gamma * max(diag(M));
            M = 1/gamma * M;
            W = eye(n) - M.M; % not subtraction-free, but the scaling ensures subtraction-safety
            Z = (1/2) * (eye(n) + M);
            p = M.v;
            for ell = 1:opts.maxit
                [L, U] = lu(Z, use_mex=opts.use_mex);
                p = p + 1/2 * W * (U\(L\p)); % p yields the triplet vector v of Z
                W = 1/4 * W * (U\(L\W)); % W is the update on the off-diagonal
                Z = mmatrix(Z.M - 1/2 * W, M.u, p + 1/2 * W * M.u);
                if norm(W) < opts.tol
                   break;
                end
            end
            out.it = ell;
            B = sqrt(gamma) * Z;
        end

    end
    %----------------------------------------------------------------------
    methods (Static)

    function [G, history] = qbd(A0, A1, A2, opts)
        % Logarithmic reduction for A0 + A1*G + A2*G^2 = G.
        %
        % G = mmatrix.qbd(A0, A1, A2)
        % [G, history] = mmatrix.qbd(A0, A1, A2, maxit=40)
        % G = mmatrix.qbd(A0, A1, A2, accurate=false)
        % G = mmatrix.qbd(A0, A1, A2, u=u, v=v)
        %
        % Inputs:
        %   A0, A1, A2: nonnegative square matrices of the same size.
        %   u, v: vectors satisfying  u > 0, v >= 0, and (I-A0-A1-A2)*u = v.
        %
        % Defaults u = ones(n,1), v = zeros(n,1) assume that A0+A1+A2 is row stochastic.
        %
        % accurate=true  uses the triplet-based GTH-like.
        % accurate=false uses MATLAB's standard linear solver.
        %
        % The method performs maxit LR updates after initialisation.
        % history{k+1} contains the iterate G_k.
        %
        % Reference:
        % Guiding Gu, Li Wang and Ren-Cang Li,
        % Highly Accurate Latouche-Ramaswami Logarithmic Reduction
        % Algorithm for Quasi-Birth-and-Death Process,
        % Journal of Mathematical Study, 55 (2022), 180-194.
        % Algorithms 3.1 and 4.1.

        arguments
            A0
            A1
            A2
            opts.u = []
            opts.v = []
            opts.maxit (1,1) double ...
                {mustBeInteger, mustBeNonnegative} = 40
            opts.accurate (1,1) logical = true
            opts.use_mex (1,1) logical = true
        end
        n = size(A0,1);
        if n == 0 || ~isequal(size(A0), [n,n]) || ~isequal(size(A1), [n,n]) || ~isequal(size(A2), [n,n])
            error('mmatrix:qbd', 'A0, A1 and A2 must be nonempty square matrices of the same size.');
        end
        if any(A0(:) < 0) || any(A1(:) < 0) || any(A2(:) < 0)
            error('mmatrix:qbd','A0, A1 and A2 must be nonnegative.');
        end
        I = eye(n);
        % Initial solve.
        if opts.accurate
            u = opts.u;
            v = opts.v;
            if isempty(u)
                u = ones(n,1);
            end
            if isempty(v)
                v = zeros(n,1);
            end
            u = u(:);
            v = v(:);
            if numel(u) ~= n || numel(v) ~= n ||any(u <= 0) || any(v < 0)
                error('mmatrix:qbd','u and v must have length n, with u > 0 and v >= 0.');
            end
            % Triplet representation of I-A1.
            M = mmatrix(-A1, u, v + A0*u + A2*u);
            [Lf, Uf] = lu(M, use_mex=opts.use_mex);
            X = Uf \ (Lf \ [A0, A2, v]);
            L = X(:,1:n);
            H = X(:,n+1:2*n);
            w = X(:,end);
            p = w + (H+L)*w;
        else
            X = (I-A1) \ [A0, A2];
            L = X(:,1:n);
            H = X(:,n+1:2*n);
        end
        G = L;
        T = H;
        history = {};
        if nargout > 1
            history = cell(opts.maxit+1,1);
            history{1} = G;
        end
        % LR updates
        for k = 1:opts.maxit
            B0 = L*L;
            B1 = H*L + L*H;
            B2 = H*H;
            if opts.accurate
                % Triplet representation of I-B1.
                M = mmatrix(-B1, u, p + B0*u + B2*u);
                [Lf, Uf] = lu(M, use_mex=opts.use_mex);
                X = Uf \ (Lf \ [B0, B2, p]);
                L = X(:,1:n);
                H = X(:,n+1:2*n);
                w = X(:,end);
                p = w + (H+L)*w;
            else
                X = (I-B1) \ [B0, B2];
                L = X(:,1:n);
                H = X(:,n+1:2*n);
            end
            G = G + T*L;
            T = T*H;
            if nargout > 1
                history{k+1} = G;
            end
        end
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
            if isnumeric(A) && isscalar(A) && A > 0 % scalar * M
                C = mmatrix(A*B.M, B.u, A*B.v, left=B.left);
            elseif isnumeric(B) && isscalar(B) && B > 0 % M * scalar
                C = mmatrix(A.M*B, A.u, A.v*B, left=A.left);
            elseif isnumeric(A) && isequal(size(A), [B.n, 1]) && all(A>0) % column * M
                % this is equivalent to diag(column) * M
                if B.left
                    C = mmatrix(A .* B.M, B.u ./ transpose(A), B.v, left=true);
                else
                    C = mmatrix(A .* B.M, B.u, A .* B.v, left=false);
                end
            elseif isnumeric(B) && isequal(size(B), [A.n, 1]) && all(B>0) % M * column
                C = times(B, A); % .* is commutative
            elseif isnumeric(A) && isequal(size(A), [1, B.n]) && all(A>0) % row * M
                % this is equivalent to M * diag(row)
                if B.left
                    C = mmatrix(A .* B.M, B.u, A .* B.v, left=true);
                else
                    C = mmatrix(A .* B.M, B.u ./ transpose(A), B.v, left=false);
                end
            elseif isnumeric(B) && isequal(size(B), [1, A.n]) && all(B>0) % M * row
                C = times(B, A); % .* is commutative
            else
                warning('mmatrix:inaccurate', 'converted to double for the product, results may not be componentwise accurate')
                C = double(A) .* double(B);
            end
        end

        function C = mtimes(A, B)
            if isscalar(A)  && A > 0
                C = mmatrix(A*B.M, B.u, A*B.v, left=B.left);
            elseif isscalar(B) && B > 0
                C = mmatrix(A.M*B, A.u, A.v*B, left=A.left);
            else
                warning('mmatrix:inaccurate', 'converted to double for the product, results may not be componentwise accurate')
                C = double(A) * double(B);
            end
        end
        function MT = transpose(M)
            MT = mmatrix(M.M', M.u', M.v', left=not(M.left));
        end
        function MT = ctranspose(M)
            MT = transpose(M);
        end
        function P = mpower(X, Y)
            if isequal(Y, 1/2)
                P = sqrtm(X);
            else
                warning('mmatrix:inaccurate', 'converted to double for the matrixpower, results may not be componentwise accurate')
                P = double(X)^double(Y);
            end
        end
        function x = norm(A, varargin)
            x = norm(double(A), varargin{:});
        end
        function d = diag(obj)
            d = transpose(obj.M(1:obj.n+1:end));
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

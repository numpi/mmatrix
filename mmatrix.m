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
            X = dmmtri(Mf);
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
            if isa(A, 'mmatrix') && isa(B, 'mmatrix') && A.left == B.left && isequal(A.u, B.u)
                C = mmatrix(A.M+B.M, A.u, A.v+B.v, left=A.left);
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
            warning('mmatrix:inaccurate', 'converted to double for the product, results may not be componentwise accurate')
            C = double(A) .* double(B);
        end
        function C = mtimes(A, B)
            warning('mmatrix:inaccurate', 'converted to double for the product, results may not be componentwise accurate')
            C = double(A) * double(B);
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

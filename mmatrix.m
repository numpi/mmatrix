classdef mmatrix
    % matrix solver functions
    
    methods(Static)
        
        %------------------------------------------------------------------
        function [A, b, w] = create_fem_problem(n, h, sigma)
            % Create FEM problem matrix and vectors
            A = full(gallery('tridiag', n, -1, 2, -1));
            A(1,1) = 1 + h*sigma;
            A(end,end) = 1 + h*sigma;
            
            w = zeros(n,1);
            w(1) = h*sigma;
            w(end) = h*sigma;
            
            b = h * sin(pi * h * (1:n)');
            b(1) = b(1)/2;
            b(end) = b(end)/2;
        end
        
        %------------------------------------------------------------------
        function x_exact = exact_solution(A, v, w, b)
            % Compute exact solution
            n = length(A);
            Avpa = vpa(A);
            Avpa(1:n+1:end) = 0;
            f = w - Avpa*v;
            Avpa(1:n+1:end) = f ./ v;
            x_exact = double(Avpa \ b);
        end
         function x_exact = exact_solution_lu(A, b)
            x_exact = double(vpa(A) \ b);
         end
        %------------------------------------------------------------------
        function [time, rel_error] = run_gthlu(A, b, x_exact, v, w)
            % Run GTH method
                tic;
                [L, U] = gthlu(A, v, w);
                time = toc;
                x = U \ (L \ b);
                rel_error = norm(x - x_exact) / norm(x_exact);
        end
        
        %------------------------------------------------------------------
        function [time, rel_error] = run_standard_lu(A, b, x_exact)
            % Run standard LU method
                tic;
                [L, U] = lu(A);
                time = toc;
                x = U \ (L \ b);
                rel_error = norm(x - x_exact) / norm(x_exact);
        end
        
        %------------------------------------------------------------------
        function [time, rel_error] = run_dmmtf2(A, b, x_exact, v, w)
            % Run DMMTF2 method
                tic;
                [Af, info] = mmex.dmmtf2(A, v, w);
                time = toc;
                    x = mmatrix.solvelu(Af, b);
                    rel_error = norm(x - x_exact) / norm(x_exact);
        end
        
        %------------------------------------------------------------------
        function [time, rel_error] = run_dmmtrf(A, b, x_exact, v, w, nb)
            % Run DMMTRF method
                tic;
                [Af, info] = mmex.dmmtrf(A, v, w, nb);
                time = toc;
                    x = mmatrix.solvelu(Af, b);
                    rel_error = norm(x - x_exact) / norm(x_exact);
        end
        
        %------------------------------------------------------------------
        function [time, rel_error] = run_dmmtrff(A, b, x_exact, v, w, nb)
            % Run DMMTRFF method
                tic;
                [Af, info] = mmex.dmmtrff(A, v, w, nb);
                time = toc;
                    x = mmatrix.solvelu(Af, b);
                    rel_error = norm(x - x_exact) / norm(x_exact);
        end
        
        %------------------------------------------------------------------
        function x = solvelu(A, b)
            % Solve using LU factors
            n = length(A);
            L = tril(A,-1) + eye(n);
            U = triu(A);
            y = L \ b;
            x = U \ y;
        end
        
        %------------------------------------------------------------------
        function results = run_dimension_experiments(params, method_flags)
            % Run experiments for varying dimensions
            sigma = params.sigma;
            ns = params.ns;
            nb = params.nb;
            
            nt = length(ns);
            results.ns = ns;
            z = zeros(1, nt);
            
            % Initialize result arrays
            results.conds = z;
            if method_flags.gth
                results.norms_gth = z;
                results.tgth = z;
            end
            if method_flags.lu
                results.norms_lu = z;
                results.tlu = z;
            end
            if method_flags.dmmtf2
                results.norms_dmmtf2 = z;
                results.tdmmtf2 = z;
            end
            if method_flags.dmmtrf
                results.norms_dmmtrf = z;
                results.tdmmtrf = z;
            end
            if method_flags.dmmtrff
                results.norms_dmmtrff = z;
                results.tdmmtrff = z;
            end
            
            for k = 1:nt
                n = ns(k);
                h = 1/(n+1);
                v = ones(n,1);
                
                % Create problem
                [A, b, w] = mmatrix.create_fem_problem(n, h, sigma);
                
                results.conds(k) = cond(A);
                
                % Exact solution
                x_exact = mmatrix.exact_solution(A, v, w, b);
                x_exact_lu = mmatrix.exact_solution_lu(A, b);

    % Properties of the class
    properties (SetAccess = public)
        off (:, :) {mustBeNonpositive} % off-diagonal of A
        u (:, 1) {mustBePositive} % positive vector	such that Au = v
        v (:, 1) {mustBeNonnegative} % nonnegative vector such that Au = v
        direction (:, :) % right = 0; left = 1
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

            [n, m] = size(off);

            if ne(n, m)
                error('offdiag(A) must be a square matrix');
            end

            if ne(size(u,1), n) || ne(size(v,1), n)
                error('u and v should have the same size as offdiag(A)');
            end

            if isequal(diag(off),zeros(n,1))
            else
                warning('A has nonzero diagonal, diagonal elements of A are not used');
                off(1:n+1:n^2) = 0;
            end

            if nargin < 4
                opts.direction = 0;
            end
            A.direction = opts.direction;

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
            if (A.direction)
                U = dmltrff_mex(A.off, A.u, A.v);
            else
                U = dmmtrff_mex(A.off, A.u, A.v);
            end
            L = tril(U, -1); L(1:n+1:end) = 1;
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
            if (A.direction)
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

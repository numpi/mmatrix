classdef MatrixSolverFunctions
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
                [Af, info] = dmmtf2_mex(A, v, w);
                time = toc;
                    x = MatrixSolverFunctions.solvelu(Af, b);
                    rel_error = norm(x - x_exact) / norm(x_exact);
        end
        
        %------------------------------------------------------------------
        function [time, rel_error] = run_dmmtrf(A, b, x_exact, v, w, nb)
            % Run DMMTRF method
                tic;
                [Af, info] = dmmtrf_mex(A, v, w, nb);
                time = toc;
                    x = MatrixSolverFunctions.solvelu(Af, b);
                    rel_error = norm(x - x_exact) / norm(x_exact);
        end
        
        %------------------------------------------------------------------
        function [time, rel_error] = run_dmmtrff(A, b, x_exact, v, w, nb)
            % Run DMMTRFF method
                tic;
                [Af, info] = dmmtrff_mex(A, v, w, nb);
                time = toc;
                    x = MatrixSolverFunctions.solvelu(Af, b);
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
                [A, b, w] = MatrixSolverFunctions.create_fem_problem(n, h, sigma);
                
                results.conds(k) = cond(A);
                
                % Exact solution
                x_exact = MatrixSolverFunctions.exact_solution(A, v, w, b);
                
                % Run selected methods
                if method_flags.gth
                    [results.tgth(k), results.norms_gth(k)] = ...
                        MatrixSolverFunctions.run_gthlu(A, b, x_exact, v, w);
                end
                
                if method_flags.lu
                    [results.tlu(k), results.norms_lu(k)] = ...
                        MatrixSolverFunctions.run_standard_lu(A, b, x_exact);
                end
                
                if method_flags.dmmtf2
                    [results.tdmmtf2(k), results.norms_dmmtf2(k)] = ...
                        MatrixSolverFunctions.run_dmmtf2(A, b, x_exact, v, w);
                end
                
                if method_flags.dmmtrf
                    [results.tdmmtrf(k), results.norms_dmmtrf(k)] = ...
                        MatrixSolverFunctions.run_dmmtrf(A, b, x_exact, v, w, nb);
                end
                
                if method_flags.dmmtrff
                    [results.tdmmtrff(k), results.norms_dmmtrff(k)] = ...
                        MatrixSolverFunctions.run_dmmtrff(A, b, x_exact, v, w, nb);
                end
            end
        end
        
        %------------------------------------------------------------------
        function display_single_results_table(results, condA)
            % Display formatted table for single matrix results
            
            fprintf('\n');
            fprintf('========================================================\n');
            fprintf('                 RESULTS \n');
            fprintf('========================================================\n');
            fprintf('Matrix Condition Number: %.2e\n\n', condA);
            
            % Create table data
            TableData = table();
            TableData.Method = results.methods';
            TableData.Error = results.errors';
            TableData.Time_s = results.times';
            
            % Display table
            disp(TableData);
        end
        
        %------------------------------------------------------------------
        function display_dimension_tables(results, params)
            % Display table 
            
            fprintf('\n');
            fprintf('================================================================================\n');
            fprintf('                   RESULTS (sigma = %.1e)\n', params.sigma);
            fprintf('================================================================================\n\n');
            
            % Create error table
            fprintf('RELATIVE ERRORS:\n');
            fprintf('%-8s', 'n');
            if isfield(results, 'norms_gth')
                fprintf(' %-12s', 'GTH');
            end
            if isfield(results, 'norms_lu')
                fprintf(' %-12s', 'LU');
            end
            if isfield(results, 'norms_dmmtf2')
                fprintf(' %-12s', 'DMMTF2');
            end
            if isfield(results, 'norms_dmmtrf')
                fprintf(' %-12s', 'DMMTRF');
            end
            if isfield(results, 'norms_dmmtrff')
                fprintf(' %-12s', 'DMMTRFF');
            end
            fprintf('\n');
            
            fprintf('%-8s', '--------');
            if isfield(results, 'norms_gth')
                fprintf(' %-12s', '------------');
            end
            if isfield(results, 'norms_lu')
                fprintf(' %-12s', '------------');
            end
            if isfield(results, 'norms_dmmtf2')
                fprintf(' %-12s', '------------');
            end
            if isfield(results, 'norms_dmmtrf')
                fprintf(' %-12s', '------------');
            end
            if isfield(results, 'norms_dmmtrff')
                fprintf(' %-12s', '------------');
            end
            fprintf('\n');
            
            for i = 1:length(results.ns)
                fprintf('%-8d', results.ns(i));
                if isfield(results, 'norms_gth')
                    fprintf(' %-12.2e', results.norms_gth(i));
                end
                if isfield(results, 'norms_lu')
                    fprintf(' %-12.2e', results.norms_lu(i));
                end
                if isfield(results, 'norms_dmmtf2')
                    fprintf(' %-12.2e', results.norms_dmmtf2(i));
                end
                if isfield(results, 'norms_dmmtrf')
                    fprintf(' %-12.2e', results.norms_dmmtrf(i));
                end
                if isfield(results, 'norms_dmmtrff')
                    fprintf(' %-12.2e', results.norms_dmmtrff(i));
                end
                fprintf('\n');
            end
            
            fprintf('\n');
            fprintf('CPU TIMES (seconds):\n');
            fprintf('%-8s', 'n');
            if isfield(results, 'tgth')
                fprintf(' %-12s', 'GTH');
            end
            if isfield(results, 'tlu')
                fprintf(' %-12s', 'LU');
            end
            if isfield(results, 'tdmmtf2')
                fprintf(' %-12s', 'DMMTF2');
            end
            if isfield(results, 'tdmmtrf')
                fprintf(' %-12s', 'DMMTRF');
            end
            if isfield(results, 'tdmmtrff')
                fprintf(' %-12s', 'DMMTRFF');
            end
            fprintf('\n');
            
            fprintf('%-8s', '--------');
            if isfield(results, 'tgth')
                fprintf(' %-12s', '------------');
            end
            if isfield(results, 'tlu')
                fprintf(' %-12s', '------------');
            end
            if isfield(results, 'tdmmtf2')
                fprintf(' %-12s', '------------');
            end
            if isfield(results, 'tdmmtrf')
                fprintf(' %-12s', '------------');
            end
            if isfield(results, 'tdmmtrff')
                fprintf(' %-12s', '------------');
            end
            fprintf('\n');
            
            for i = 1:length(results.ns)
                fprintf('%-8d', results.ns(i));
                if isfield(results, 'tgth')
                    fprintf(' %-12.4f', results.tgth(i));
                end
                if isfield(results, 'tlu')
                    fprintf(' %-12.4f', results.tlu(i));
                end
                if isfield(results, 'tdmmtf2')
                    fprintf(' %-12.4f', results.tdmmtf2(i));
                end
                if isfield(results, 'tdmmtrf')
                    fprintf(' %-12.4f', results.tdmmtrf(i));
                end
                if isfield(results, 'tdmmtrff')
                    fprintf(' %-12.4f', results.tdmmtrff(i));
                end
                fprintf('\n');
            end
            
            fprintf('\n');
            fprintf('CONDITION NUMBERS:\n');
            fprintf('%-8s', 'n');
            for i = 1:length(results.ns)
                fprintf(' %-12d', results.ns(i));
            end
            fprintf('\n');
            fprintf('%-8s', 'cond(A)');
            for i = 1:length(results.ns)
                fprintf(' %-12.2e', results.conds(i));
            end
            fprintf('\n');
        end
    end
end
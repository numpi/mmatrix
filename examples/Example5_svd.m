% Example 5: Accuracy of SVD using a finite-element matrix
%
%   This script compares standard MATLAB SVD and the mmatrix-toolbox SVD.
%
%  Notes:
%   * Standard MATLAB SVD loses component-wise accuracy (compS_A) as the 
%     matrix condition number increases.
%   * The SVD algorithm in the mmatrix-toolbox maintains componentwise 
%     errors (compS_M) close to machine precision.
%   * The norm-wise error is dominated by the largest singular values. 
%     As a result, it stays near machine precision for both methods and 
%     hides significant errors occurring in the smaller singular values.
clear; clc;

sigmas = [1e-1, 1e-6, 1e-10];  
matrix_sizes = [5, 16, 32]; % Matrix dimensions (n x n)
num_sizes = length(matrix_sizes);

for sigma = sigmas    
    condA    = zeros(num_sizes, 1); % Condition number of the FEM matrix
    compS_A  = zeros(num_sizes, 1); % Componentwise error (Standard MATLAB SVD)
    normS_A  = zeros(num_sizes, 1); % Normwise error (Standard MATLAB SVD)
    compS_M = zeros(num_sizes, 1);  % Componentwise error (mmatrix-toolbox SVD)
    normS_M = zeros(num_sizes, 1);  % Normwise error (mmatrix-toolbox SVD)

    for k = 1:num_sizes
        n = matrix_sizes(k);
        h = 1 / (n + 1);
        
        % Generate the FEM Matrix
        [M, ~] = gettriplet('fem', size=n, left=0, sigma=sigma, gamma=1, h=h);
        A = M.M;         
        condA(k) = cond(A);
        
        % High-Precision SVD
        [~, S_h, ~] = svd(vpa(M, 100));
        S_h = (diag(S_h)); 

        % Standard MATLAB SVD
        [~, S_m, ~] = svd(A);  
        S_m = diag(S_m);
        
        % Compute errors for standard SVD
        compS_A(k) = max(abs((S_m - S_h) ./ S_h));
        normS_A(k) = norm(S_m - S_h) / norm(S_h);
        
        % mmatrix-toolbox SVD
        [~, S, ~] = svd(M);  
        S = diag(S);
        
        % Compute errors for mmatrix-toolbox SVD
        compS_M(k) = max(abs((S - S_h) ./ S_h));
        normS_M(k) = norm(S - S_h) / norm(S_h);
        
        fprintf('Matrix size n = %2d processed (cond(A) = %.2e)\n', n, condA(k));
    end
    
    ResultsTable = table(matrix_sizes', condA, compS_A, normS_A, compS_M, normS_M, ...
        'VariableNames', {'n', 'condA', 'compS_A', 'normS_A', 'compS_M', 'normS_M'});
    
    disp(ResultsTable);
end

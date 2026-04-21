function [M, out] = gettriplet(m, opts)
% Load examples of M-matrices with triplet
% M = gettriplet(m)
% [M, out] = gettriplet(m, opts)
% 
% n: example
% in parentheses the default value of the options
% opts.size: size of the data (10)
% opts.left: left triplet (0)
% M: the example
% out: out structure with information
% 
% 'rand': M-matrix with uniformly distributed component, off-diagonal in
% [-1,0], u in [0,1], v in [0, opts.sigma (1)], 
% rng initialized with opts.seed ("shuffle")
%
% 'fem': example from Finite Element Methods with step opts.h

arguments
    m
    opts.size = 10
    opts.left = 0
    opts.sigma = 1
    opts.gamma = 1
    opts.h = 1d-2;
    opts.seed = 'shuffle';
end
n = opts.size;

switch m
    
    case {'rand', 1}
        sigma = opts.sigma;
        seed = opts.seed;
        rng(seed); % rng selection for reproducibility
        M = mmatrix(-rand(n), rand(n,1), sigma * rand(n,1), ...
                    left = opts.left);
    
    case {'fem', 2}
        if opts.left
            warning('This example has no left triplet, right is used')
        end
        sigma = opts.sigma;
        h = opts.h;
        M = full(gallery('tridiag', n, -1, 2, -1));
        M(1, 1) = 1 + h * sigma;
        M(n, n) = 1 + h * sigma;
        u = ones(n, 1);
        v = zeros(n, 1);
        v([1 n]) = h * sigma * [1 1];
        M = mmatrix(M, u, v);

        b = h * sin(pi*h*(1:n)');
        b([1 n]) = b([1 n])/2;
        out.b = b;
    
    case 'dario3'
        M = - triu(ones(n));
        M(1:n+1:n^2) = n;
        for ell = 2:floor(n/4)
            M(ell:n+1:n^2) = -1;
        end
        u = ones(n, 1);
        v = [zeros(floor(n/4),1); (1:n-floor(n/4))'];
        M = mmatrix(M, u, v);

end

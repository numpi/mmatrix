function [M, out] = gettriplet(m, opts)
% Load examples of M-matrices with triplet
% M = gettriplet(m)
% [M, out] = gettriplet(m, opts)
% 
% n: example
% opts.size: size of the data (default = 10)
% opts.left: left or right triplet (default right)
% M: the example
% out: out structure with information
% 
% 'rand': M-matrix with uniformly distributed component, off-diagonal in
% [-1,0], u in [0,1], v in [0, opts.sigma]
%
% 'fem': example from Finite Element Methods with step opts.h

arguments
    m
    opts.size = 10
    opts.left = 0
    opts.sigma = 1
    opts.gamma = 1
    opts.h = 1d-2;
end
n = opts.size;

switch m
    case {'rand', 1}
        sigma = opts.sigma;
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
end

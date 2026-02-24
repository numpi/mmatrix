function [A, v, w] = randtriplet(n, opts)
% construct random m-matrix (convenient for testing)
arguments
    n
    opts.sigma = 1;
    opts.left = false;
end

A = -rand(n, n);
v = rand(n, 1);
w = opts.sigma * rand(n, 1);
% adjust diagonal
A(1:n+1:end) = 0;
A(1:n+1:end) = (w - A*v) ./ v;
if opts.left
    A = A';
    v = v';
    w = w';
end
end

function tests = BasicSolverTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Add ../src to path for the duration of this test file
import matlab.unittest.fixtures.PathFixture
testCase.applyFixture(PathFixture('..'));
end

function testGthLU(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    [L, U, ww] = gthlu(A, v, w);
    testCase.verifyEqual(L*U, A, RelTol=1e-9);
    testCase.verifyEqual(L, tril(L));
    testCase.verifyEqual(diag(L), ones(n,1));
    testCase.verifyEqual(U, triu(U));    
    testCase.verifyEqual(ww, L \ w);    

    [A, v, w] = randtriplet(n, left=true);
    [L, U, ww] = gthlu(A, v, w, 'left');
    testCase.verifyEqual(L*U, A, RelTol=1e-9);
    testCase.verifyEqual(diag(L), ones(n,1));
    testCase.verifyEqual(ww, w / U);    
end

function testDmmtf2(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    LU = dmmtf2_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    % test with more outputs
    [LU, info] = dmmtf2_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);    
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    testCase.verifyEqual(info, 0);
    [LU, info, outw] = dmmtf2_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    testCase.verifyEqual(info, 0);
    testCase.verifyEqual(L\w, outw, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = dmmtf2_mexsolve(A, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);

    % test all outputs from _mexsolve
    [x, LU, info, outw] = dmmtf2_mexsolve(A, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    testCase.verifyEqual(info, 0);
    testCase.verifyEqual(L\w, outw, 'RelTol', 1e-9);
end

function testDmmtrf2(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    LU = dmmtrf2_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = dmmtrf2_mexsolve(A, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmmtrf(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    % try different block sizes
    for blocksize = [1 2 n/2 n-1 n n+1 2*n]
        LU = dmmtrf_mex(A, v, w, blocksize);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    end
    % test also with the default block size
    LU = dmmtrf_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = dmmtrf_mexsolve(A, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmmtrff(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    % try different block sizes
    for blocksize = [1 2 n/2 n-1 n n+1 2*n]
        LU = dmmtrff_mex(A, v, w, blocksize);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    end
    % test also with the default block size
    LU = dmmtrff_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = dmmtrff_mexsolve(A, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDgetf2(testCase)
    rng(0);
    n = 4;
    [A, ~, ~] = randtriplet(n);
    [LU, ipiv] = dgetf2_mex(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);

    b = rand(n, 2);
    x = dgetf2_mexsolve(A, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDgetrf2(testCase)
    rng(0);
    n = 4;
    [A, ~, ~] = randtriplet(n);
    [LU, ipiv] = dgetrf2_mex(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);

    b = rand(n, 2);
    x = dgetrf2_mexsolve(A, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDgetrf(testCase)
    rng(0);
    n = 4;
    [A, ~, ~] = randtriplet(n);
    [LU, ipiv] = dgetrf_mex(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);

    b = rand(n, 2);
    x = dgetrf_mexsolve(A, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmltf2(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n, left=true);
    LU = dmltf2_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    % test with more outputs
    [LU, info] = dmltf2_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);    
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    testCase.verifyEqual(info, 0);
    [LU, info, outw] = dmltf2_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    testCase.verifyEqual(info, 0);
    testCase.verifyEqual(w / U, outw, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = dmltf2_mexsolve(A, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmltrf2(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n, left=true);
    [LU, info, outw] = dmltrf2_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    testCase.verifyEqual(info, 0);
    testCase.verifyEqual(w / U, outw, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = dmltrf2_mexsolve(A, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmltrf(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n, left=true);
    % try different block sizes
    for blocksize = [1 2 n/2 n-1 n n+1 2*n]
        [LU, info, outw] = dmltrf_mex(A, v, w, blocksize);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        testCase.verifyEqual(info, 0);
        % testCase.verifyEqual(w / U, outw, 'RelTol', 1e-9);
    end
    % test also with the default block size
    LU = dmltrf_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = dmltrf_mexsolve(A, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmltrff(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n, left=true);
    % try different block sizes
    for blocksize = [1 2 n/2 n-1 n n+1 2*n]
        [LU, info, outw] = dmltrff_mex(A, v, w, blocksize);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        testCase.verifyEqual(info, 0);
        % testCase.verifyEqual(w / U, outw, 'RelTol', 1e-9);
    end
    % test also with the default block size
    LU = dmltrff_mex(A, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = dmltrff_mexsolve(A, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);    
end


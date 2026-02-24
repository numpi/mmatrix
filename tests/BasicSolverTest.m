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
    [A, v, w] = randtriplet(4);
    [L, U] = gthlu(A, v, w);
    testCase.verifyEqual(L*U, A, RelTol=1e-9);
    testCase.verifyEqual(L, tril(L));
    testCase.verifyEqual(U, triu(U));    
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
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    [LU, info, outw] = dmmtf2_mex(A, v, w);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
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
end

function testDgetf2(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    [LU, ipiv] = dgetf2_mex(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);
end

function testDgetrf2(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    [LU, ipiv] = dgetrf2_mex(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);
end

function testDgetrf(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    [LU, ipiv] = dgetrf_mex(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);
end

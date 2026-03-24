function tests = BasicSolverTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Add .. to path for the duration of this test file
import matlab.unittest.fixtures.PathFixture
testCase.applyFixture(PathFixture('..'));
end

function testGthLU(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    A0 = A; A0(1:n+1:end) = nan; % sets diagonal to NaN to make sure we are not cheating
    [L, U, ww] = gthlu(A0, v, w);
    testCase.verifyEqual(L*U, A, RelTol=1e-9);
    testCase.verifyEqual(L, tril(L));
    testCase.verifyEqual(diag(L), ones(n,1));
    testCase.verifyEqual(U, triu(U));    
    testCase.verifyEqual(ww, L \ w);    

    [A, v, w] = randtriplet(n, left=true);
    A0 = A; A0(1:n+1:end) = 0;
    [L, U, ww] = gthlu(A0, v, w, 'left');
    testCase.verifyEqual(L*U, A, RelTol=1e-9);
    testCase.verifyEqual(diag(L), ones(n,1));
    testCase.verifyEqual(ww, w / U);    
end

function testDmmtf2(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    A0 = A; A0(1:n+1:end) = nan;
    LU = mmex.dmmtf2(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    % test with more outputs
    [LU, info] = mmex.dmmtf2(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);    
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    testCase.verifyEqual(info, 0);
    [LU, info, outw] = mmex.dmmtf2(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    testCase.verifyEqual(info, 0);
    testCase.verifyEqual(L\w, outw, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dmmtf2_solve(A0, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);

    % test all outputs from _solve
    [x, LU, info, outw] = mmex.dmmtf2_solve(A0, v, w, b);
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
    A0 = A; A0(1:n+1:end) = nan;
    LU = mmex.dmmtrf2(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dmmtrf2_solve(A0, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmmtrf(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    A0 = A; A0(1:n+1:end) = nan;
    % try different block sizes
    for blocksize = [1 2 n/2 n-1 n n+1 2*n]
        LU = mmex.dmmtrf(A0, v, w, blocksize);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    end
    % test also with the default block size
    LU = mmex.dmmtrf(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dmmtrf_solve(A0, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmmtrff(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n);
    A0 = A; A0(1:n+1:end) = nan;
    % try different block sizes
    for blocksize = [1 2 n/2 n-1 n n+1 2*n]
        LU = mmex.dmmtrff(A0, v, w, blocksize);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    end
    % test also with the default block size
    LU = mmex.dmmtrff(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dmmtrff_solve(A0, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDgetf2(testCase)
    rng(0);
    n = 4;
    [A, ~, ~] = randtriplet(n);
    [LU, ipiv] = mmex.dgetf2(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dgetf2_solve(A, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDgetrf2(testCase)
    rng(0);
    n = 4;
    [A, ~, ~] = randtriplet(n);
    [LU, ipiv] = mmex.dgetrf2(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dgetrf2_solve(A, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDgetrf(testCase)
    rng(0);
    n = 4;
    [A, ~, ~] = randtriplet(n);
    [LU, ipiv] = mmex.dgetrf(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dgetrf_solve(A, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDgetf2_reference(testCase)
    rng(0);
    n = 4;
    [A, ~, ~] = randtriplet(n);
    [LU, ipiv] = mmex.dgetf2_reference(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dgetf2_reference_solve(A, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDgetrf2_reference(testCase)
    rng(0);
    n = 4;
    [A, ~, ~] = randtriplet(n);
    [LU, ipiv] = mmex.dgetrf2_reference(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dgetrf2_reference_solve(A, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDgetrf_reference(testCase)
    rng(0);
    n = 4;
    [A, ~, ~] = randtriplet(n);
    [LU, ipiv] = mmex.dgetrf_reference(A);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dgetrf_reference_solve(A, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmltf2(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n, left=true);
    A0 = A; A0(1:n+1:end) = nan;
    LU = mmex.dmltf2(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    % test with more outputs
    [LU, info] = mmex.dmltf2(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);    
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    testCase.verifyEqual(info, 0);
    [LU, info, outw] = mmex.dmltf2(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    testCase.verifyEqual(info, 0);
    testCase.verifyEqual(w / U, outw, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dmltf2_solve(A0, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmltrf2(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n, left=true);
    A0 = A; A0(1:n+1:end) = nan;
    [LU, info, outw] = mmex.dmltrf2(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
    testCase.verifyEqual(info, 0);
    testCase.verifyEqual(w / U, outw, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dmltrf2_solve(A0, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmltrf(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n, left=true);
    A0 = A; A0(1:n+1:end) = nan;
    % try different block sizes
    for blocksize = [1 2 n/2 n-1 n n+1 2*n]
        [LU, info, outw] = mmex.dmltrf(A0, v, w, blocksize);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        testCase.verifyEqual(info, 0);
        % testCase.verifyEqual(w / U, outw, 'RelTol', 1e-9);
    end
    % test also with the default block size
    LU = mmex.dmltrf(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dmltrf_solve(A0, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
end

function testDmltrff(testCase)
    rng(0);
    n = 4;
    [A, v, w] = randtriplet(n, left=true);
    A0 = A; A0(1:n+1:end) = nan;
    % try different block sizes
    for blocksize = [1 2 n/2 n-1 n n+1 2*n]
        [LU, info, outw] = mmex.dmltrff(A0, v, w, blocksize);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        testCase.verifyEqual(info, 0);
        % testCase.verifyEqual(w / U, outw, 'RelTol', 1e-9);
    end
    % test also with the default block size
    LU = mmex.dmltrff(A0, v, w);
    L = eye(n) + tril(LU, -1);
    U = triu(LU);
    testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);

    b = rand(n, 2);
    x = mmex.dmltrff_solve(A0, v, w, b);
    testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);    
end


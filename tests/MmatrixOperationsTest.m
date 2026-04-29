function tests = MmatrixOperationsTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Add .. to path for the duration of this test file
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('..'));
end

function testCreation(testCase)
    for left = [true, false]
        rng(0);
        n = 7;
        AM = -rand(n, n);
        u = rand(n, 1);
        v = rand(n, 1);
        if left
            u = u'; v = v';
        end
        A = mmatrix(AM, u, v, left=left);
        if left
            testCase.verifyEqual(u*double(A), v, RelTol=1e-9);
        else
            testCase.verifyEqual(double(A)*u, v, RelTol=1e-9);
        end

        % test constructing without v
        testCase.verifyWarning(@() mmatrix(double(A), u, left=left), 'mmatrix:inaccurate');
        DM = diag(diag(A));
        testCase.verifyWarningFree(@() mmatrix(DM, u, left=left), 'mmatrix:inaccurate');
        D = mmatrix(DM, u, left=left);
        testCase.verifyEqual(DM, double(D), RelTol=1e-9);

        state = warning;
        warning('off', 'mmatrix:inaccurate');
            AA = mmatrix(double(A), u, left=left);
            testCase.verifyEqual(double(A), double(AA), RelTol=1e-9);
        warning(state);
    end
end

function testWarnings(testCase)
    for left = [true, false]
        rng(0);
        n = 5;
        A = gettriplet('rand', size=n, left=left, seed=1);
        B = gettriplet('rand', size=n, left=false, seed=2);
        C = gettriplet('rand', size=n, left=true, seed=3);

        D = rand(n);
        state = warning;
        warning('on', 'mmatrix:inaccurate');
        % mmatrix op double
        testCase.verifyWarning(@() A+D, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() A-D, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() A*D, 'mmatrix:inaccurate');
        % double op mmatrix
        testCase.verifyWarning(@() D+A, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() D-A, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() D*A, 'mmatrix:inaccurate');
        % mmatrix op mmatrix
        testCase.verifyWarning(@() A+B, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() A-B, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() A*B, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() A+C, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() A-C, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() A*C, 'mmatrix:inaccurate');
        % systems with indefinite rhs
        D(1:2,1) = [1 -1];
        testCase.verifyWarning(@() A \ D, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() D / A, 'mmatrix:inaccurate');
        warning(state);
    end
end

function testSum(testCase)
    for left = [true, false]
        rng(0);
        n = 5;
        [AM, u, v] = randtriplet(n, left=left);
        [BM, uB, vB] = randtriplet(n, left=left);
    
        A = mmatrix(AM, u, v, left=left);
        B = mmatrix(BM, uB, vB, left=left);
    
        testCase.verifyWarning(@() A+B, 'mmatrix:inaccurate')
        state = warning;
        warning('off', 'mmatrix:inaccurate');
            C = A + B;
            testCase.verifyClass(C, 'double');
            testCase.verifyEqual(C, double(A) + double(B), RelTol=1e-9);
        warning(state);

        % modify B so that it has uB==uA
        BM = BM .* transpose((uB./u));
        if left
            testCase.verifyEqual(u*BM, vB, RelTol=1e-9);
        else
            testCase.verifyEqual(BM*u, vB, RelTol=1e-9);
        end
        uB = u;
        B = mmatrix(BM, uB, vB, left=left);
    
        C = A + B;
        testCase.verifyClass(C, 'mmatrix');
        testCase.verifyEqual(double(C), double(A)+double(B), RelTol=1e-9);
    end
end

function testmldivide(testCase)
    for left = [true, false]
        rng(0);
        n = 9;
        [AM, u, v] = randtriplet(n, left=left);
        A = mmatrix(AM, u, v, left=left);
        D = rand(n, 4);
        testCase.verifyEqual(A \ D, double(A) \ D, RelTol=1e-9);
        testCase.verifyEqual(D' / A, D' / double(A), RelTol=1e-9);
    end
end

function testUtilityFunctions(testCase)
    for left = [true, false]
        rng(0);
        n = 9;
        [AM, u, v] = randtriplet(n, left=left);
        A = mmatrix(AM, u, v, left=left);
        testCase.verifyEqual(size(A), [n n]);
        testCase.verifyEqual(size(A,1), n);
        testCase.verifyEqual(length(A), n);
        testCase.verifyEqual(norm(A), norm(double(A)), RelTol=1e-9);
        testCase.verifyEqual(norm(A, 1), norm(double(A), 1), RelTol=1e-9);
    end
end

function testIndexing(testCase)
    for left = [true, false]
        rng(0);
        n = 3;
        [AM, u, v] = randtriplet(n, left=left);
        A = mmatrix(AM, u, v, left=left);
        % Test indexing for mmatrix
        I1 = [false true true];
        I2 = 1:2;
        testCase.verifyEqual(A(I1, I2), AM(I1,I2), RelTol=1e-9);
        testCase.verifyEqual(double(A(I1, I1)), AM(I1,I1), RelTol=1e-9);
        testCase.verifyClass(A(I1, I1), 'mmatrix');
        testCase.verifyEqual(double(A(I2, I2)), AM(I2,I2), RelTol=1e-9);
        testCase.verifyClass(A(I2, I2), 'mmatrix');
    end
end

function testTimes(testCase)
    for left = [true, false]
        n = 5;
        A = gettriplet('rand', size=n, seed=0, left=left);
        alpha = 4.5;
        v = rand(n, 1);
        testCase.verifyEqual(double(A .* alpha), double(A) .* alpha, RelTol=1e-9);
        testCase.verifyClass(A .* alpha, 'mmatrix');
        testCase.verifyEqual(double(alpha .* A), alpha .* double(A), RelTol=1e-9);
        testCase.verifyClass(alpha .* A, 'mmatrix');
        testCase.verifyEqual(double(A .* v), double(A) .* v, RelTol=1e-9);
        testCase.verifyClass(A .* v, 'mmatrix');
        testCase.verifyEqual(double(v .* A), v .* double(A), RelTol=1e-9);
        testCase.verifyClass(v .* A, 'mmatrix');
        testCase.verifyEqual(double(A .* v'), double(A) .* v', RelTol=1e-9);
        testCase.verifyClass(A .* v', 'mmatrix');
        testCase.verifyEqual(double(v' .* A), v' .* double(A), RelTol=1e-9);
        testCase.verifyClass(v' .* A, 'mmatrix');

        % some cases that throw warnings...
        testCase.verifyWarning(@() (-alpha) .* A, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() (-v) .* A, 'mmatrix:inaccurate');
        testCase.verifyWarning(@() (-v') .* A, 'mmatrix:inaccurate');
        % ...and return double
        state = warning;
        warning('off', 'mmatrix:inaccurate');
        testCase.verifyClass((-alpha) .* A, 'double');
        testCase.verifyClass((-v) .* A, 'double');
        testCase.verifyClass((-v') .* A, 'double');
        warning(state);
    end
end

function testLU(testCase)
    for left = [true, false]
        rng(0);
        n = 8;
        [AM, u, v] = randtriplet(n, left=left);
        A = mmatrix(AM, u, v, left=left);
        % Test indexing for mmatrix
        [L, U] = lu(A);
        testCase.verifyEqual(double(A), L*U, RelTol=1e-9);        
        [L, U] = lumatlab(A);
        testCase.verifyEqual(double(A), L*U, RelTol=1e-9);        
    end
end

function testTranspose(testCase)
    for left = [true, false]
        M = gettriplet('rand', seed=0, left=left);
        testCase.verifyEqual(double(transpose(M)), transpose(double(M)), RelTol=1e-9);
        testCase.verifyEqual(double(M'), (double(M))', RelTol=1e-9);
    end
end

function testDiag(testCase)
    M = gettriplet('rand', seed=0);
    testCase.verifyEqual(diag(double(M)), diag(M));
end
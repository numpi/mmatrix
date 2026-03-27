function tests = BasicSolverTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Add .. to path for the duration of this test file
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('..'));
end

function testSum(testCase)
    rng(0);
    n = 5;
    [AM, u, v] = randtriplet(n);
    [BM, uB, vB] = randtriplet(n);

    A = mmatrix(AM, u, v);
    B = mmatrix(BM, uB, vB);

    state = warning;
    testCase.verifyWarning(@() A+B, 'mmatrix:inaccurate')
    warning('off', 'mmatrix:inaccurate');
        C = A + B;
        testCase.verifyClass(C, 'double');
        testCase.verifyEqual(C, double(A) + double(B), RelTol=1e-9);
    warning(state);

    % modify B so that it has uB==uA
    BM = BM .* transpose((uB./u));
    uB = u;
    testCase.verifyEqual(BM*uB, vB, RelTol=1e-9);
    B = mmatrix(BM, uB, vB);

    C = A + B;
    testCase.verifyClass(C, 'mmatrix');
    testCase.verifyEqual(double(C), double(A)+double(B), RelTol=1e-9);
end

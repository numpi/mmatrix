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
        [AM, u, v] = randtriplet(n, left=left);
        A = mmatrix(AM, u, v, left=left);
        testCase.verifyEqual(AM, double(A), RelTol=1e-9);
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
    
        state = warning;
        testCase.verifyWarning(@() A+B, 'mmatrix:inaccurate')
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

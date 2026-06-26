function tests = EigMinTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Add .. to path for the duration of this test file
import matlab.unittest.fixtures.PathFixture
testCase.applyFixture(PathFixture('..'));
end

function EigMinBasicTest(testCase)
    for left = [true, false]
        n = 7;
        M = gettriplet('rand', size=n, seed=0, sigma=0.2, left=left);
        testCase.verifyEqual(eigmin(M), min(eig(double(M))), RelTol=1e-9);
    end
end

function EigMinReducibleZeroTest(testCase)
    A = mmatrix([1 -1; 0 0], [1;1], [0;0]);
    testCase.verifyEqual(eigmin(A), 0);
end
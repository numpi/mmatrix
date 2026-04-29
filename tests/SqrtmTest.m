function tests = SqrtmTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Add .. to path for the duration of this test file
import matlab.unittest.fixtures.PathFixture
testCase.applyFixture(PathFixture('..'));
end

function SqrtmBasicTest(testCase)
    for left = [true, false]
        n = 7;
        M = gettriplet('rand', size=n, seed=0, sigma=0.2, left=left);
        S = sqrtm(M);
        testCase.verifyEqual(double(S)^2, double(M), RelTol=1e-9);
    end
end

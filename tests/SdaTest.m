function tests = SDATest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Add .. to path for the duration of this test file
import matlab.unittest.fixtures.PathFixture
testCase.applyFixture(PathFixture('..'));
end

function SDABasicTest(testCase)
    for left = [true, false]
        n = 7;
        n1 = 3;
        M = gettriplet('rand', size=n, seed=0, sigma=0.2, left=left);
        B = M(1:n1, 1:n1);
        D = -M(1:n1, n1+1:end);
        C = -M(n1+1:end, 1:n1);
        A = M(n1+1:end, n1+1:end);
        [X, Y] = sda(M, n1);
    
        testCase.verifyEqual(X*D*X + C, double(A)*X + X*double(B), RelTol=1e-9);
        testCase.verifyEqual(D + Y*C*Y, Y*double(A) + double(B)*Y, RelTol=1e-9);
    end
end

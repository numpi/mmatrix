function tests = RandTests
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Add ../src to path for the duration of this test file
import matlab.unittest.fixtures.PathFixture
testCase.applyFixture(PathFixture('..'));
end

function testRandtriplet(testCase)
    n = 4;
    [A, v, w] = randtriplet(n);
    offdiag = logical(ones(n) - eye(n));
    testCase.verifyLessThan(A(offdiag), 0);
    testCase.verifyGreaterThan(v, 0);
    testCase.verifyGreaterThan(w, 0);
    testCase.verifyEqual(A*v, w, RelTol=1e-9);
end

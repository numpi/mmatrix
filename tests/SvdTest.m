function tests = SvdTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Add .. to path for the duration of this test file
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture('..'));
end

function testBasicSVD(testCase)
    for left = [true, false]
        n = 7;
        M = gettriplet('rand', size=n, left=left);
        [U, S, V] = svd(M);
        testCase.verifyEqual(U'*U, eye(n), AbsTol=1e-9);
        testCase.verifyEqual(V'*V, eye(n), AbsTol=1e-9);
        testCase.verifyEqual(U*S*V', double(M), RelTol=1e-9);
        testCase.verifyEqual(svd(M), diag(S));
    end    
end

% TODO: test componentwise accuracy in a critical case

function tests = SchurTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% Add .. to path for the duration of this test file
import matlab.unittest.fixtures.PathFixture
testCase.applyFixture(PathFixture('..'));
end

function testSchur(testCase)
    n = 7;
    A = gettriplet('rand', size=n, seed=0);
    I = [true, false, true, false, false, true, true];
    notI = not(I);
    J = 3:7;
    notJ = 1:2;
    DA = double(A);
    testCase.verifyEqual(double(schur(A, I)), DA(notI,notI) - DA(notI,I)/DA(I,I)*DA(I,notI), RelTol=1e-9);
    testCase.verifyEqual(double(schur(A, J)), DA(notJ,notJ) - DA(notJ,J)/DA(J,J)*DA(J,notJ), RelTol=1e-9);
end
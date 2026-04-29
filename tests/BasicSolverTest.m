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
    M = gettriplet('rand', size=n);
    A = M.M; v = M.u; w = M.v;
    A0 = A; A0(1:n+1:end) = nan; % sets diagonal to NaN to make sure we are not cheating
    [L, U, ww] = gthlu(A0, v, w);
    testCase.verifyEqual(L*U, A, RelTol=1e-9);
    testCase.verifyEqual(L, tril(L));
    testCase.verifyEqual(diag(L), ones(n,1));
    testCase.verifyEqual(U, triu(U));    
    testCase.verifyEqual(ww, L \ w);    

    M = gettriplet('rand', size=n, left=true);
    A = M.M; v = M.u; w = M.v;
    A0 = A; A0(1:n+1:end) = 0;
    [L, U, ww] = gthlu(A0, v, w, 'left');
    testCase.verifyEqual(L*U, A, RelTol=1e-9);
    testCase.verifyEqual(diag(L), ones(n,1));
    testCase.verifyEqual(ww, w / U);
end

function testRightTripletSolvers(testCase)
    rng(0);
    n = 4;
    M = gettriplet('rand', size=n);
    A = M.M; v = M.u; w = M.v;
    A0 = A; A0(1:n+1:end) = nan;
    solvers = { "dmmtf2", "dmmtrf2", "dmmtrf", "dmmtrff" };
    for k = 1:length(solvers)
        solver = str2func("mmex."+solvers{k});
        LU = solver(A0, v, w);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        % test with more outputs
        [LU, info] = solver(A0, v, w);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);    
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        testCase.verifyEqual(info, 0);
        [LU, info, outw] = solver(A0, v, w);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        testCase.verifyEqual(info, 0);
        testCase.verifyEqual(L\w, outw, 'RelTol', 1e-9);
        % test _solve        
        solver = str2func("mmex."+solvers{k} + "_solve");
        b = rand(n, 2);
        x = solver(A0, v, w, b);
        testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
        % test all outputs from _solve
        [x, LU, info, outw] = solver(A0, v, w, b);
        testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        testCase.verifyEqual(info, 0);
        testCase.verifyEqual(L\w, outw, 'RelTol', 1e-9);
        % test _lsolve
        solver = str2func("mmex."+solvers{k} + "_lsolve");
        b = rand(n, 2);
        x = solver(A0, v, w, b');
        testCase.verifyEqual(x, b'/A, 'RelTol', 1e-9);
    end
end    

function testLeftTripletSolvers(testCase)
    rng(0);
    n = 4;
    M = gettriplet('rand', size=n, left=true);
    A = M.M; v = M.u; w = M.v;
    A0 = A; A0(1:n+1:end) = nan;
    solvers = {"dmltf2", "dmltrf2", "dmltrf", "dmltrff"};
    for k = 1:length(solvers)
        solver = str2func("mmex."+solvers{k});
        LU = solver(A0, v, w);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        % test with more outputs
        [LU, info] = solver(A0, v, w);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);    
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        testCase.verifyEqual(info, 0);
        [LU, info, outw] = solver(A0, v, w);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        testCase.verifyEqual(info, 0);
        testCase.verifyEqual(w / U, outw, 'RelTol', 1e-9);
        % test _solve        
        solver = str2func("mmex."+solvers{k}+"_solve");
        b = rand(n, 2);
        x = solver(A0, v, w, b);
        testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
        % test all outputs from _solve
        [x, LU, info, outw] = solver(A0, v, w, b);
        testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
        testCase.verifyEqual(info, 0);
        testCase.verifyEqual(w / U, outw, 'RelTol', 1e-9);
        % test _lsolve
        solver = str2func("mmex."+solvers{k} + "_lsolve");
        b = rand(n, 2);
        x = solver(A0, v, w, b');
        testCase.verifyEqual(x, b'/A, 'RelTol', 1e-9);
    end
end

function testTraditionalSolvers(testCase)
    rng(0);
    n = 4;
    M = gettriplet('rand', size=n, seed=0);
    A = M.M; % v = M.u; w = M.v;
    solvers = {"dgetf2", "dgetrf2", "dgetrf", "dgetf2_reference", "dgetrf2_reference", "dgetrf_reference"};
    for k = 1:length(solvers)
        solver = str2func("mmex."+solvers{k});
        [LU, ipiv] = solver(A);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);
        % test _solve        
        solver = str2func("mmex."+solvers{k}+"_solve");
        b = rand(n, 2);
        x = solver(A, b);
        testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
        % test all outputs from _solve
        [x, LU, ipiv, info] = solver(A, b);
        testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
        L = eye(n) + tril(LU, -1);
        U = triu(LU);
        testCase.verifyEqual(L*U, A(ipiv, :), 'RelTol', 1e-9);
        testCase.verifyEqual(info, 0);
        % Commented because _lsolve is not implemented for these
        % test _lsolve
        % solver = str2func("mmex."+solvers{k} + "_lsolve");
        % b = rand(n, 2);
        % x = solver(A, b');
        % testCase.verifyEqual(x, b'/A, 'RelTol', 1e-9);
    end
end

function testRightTripletSolversWithBlockSizes(testCase)
    rng(0);
    n = 6;
    M = gettriplet('rand', size=n);
    A = M.M; v = M.u; w = M.v;
    A0 = A; A0(1:n+1:end) = nan;
    b = rand(n, 3);
    solvers = {"dmmtrf", "dmmtrff" };
    for k = 1:length(solvers)        
        for blocksize = [1 2 n/2 n-1 n n+1 2*n]
            solver = str2func("mmex."+solvers{k});
            [LU, info, outw] = solver(A0, v, w, blocksize);
            L = eye(n) + tril(LU, -1);
            U = triu(LU);
            testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);            
            testCase.verifyEqual(info, 0);
            testCase.verifyEqual(L\w, outw, 'RelTol', 1e-9);
            solver = str2func("mmex."+solvers{k} + "_solve");            
            x = solver(A0, v, w, b, blocksize);
            testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
            solver = str2func("mmex."+solvers{k} + "_lsolve");
            x = solver(A0, v, w, b', blocksize);
            testCase.verifyEqual(x, b'/A, 'RelTol', 1e-9);
        end    
    end
end

function testLeftTripletSolversWithBlockSizes(testCase)
    rng(0);
    n = 6;
    M = gettriplet('rand', size=n, left=true);
    A = M.M; v = M.u; w = M.v;
    A0 = A; A0(1:n+1:end) = nan;
    b = rand(n, 3);
    solvers = {"dmltrf", "dmltrff"};
    for k = 1:length(solvers)        
        for blocksize = [1 2 n/2 n-1 n n+1 2*n]
            solver = str2func("mmex."+solvers{k});            
            [LU, info, outw] = solver(A0, v, w, blocksize);
            L = eye(n) + tril(LU, -1);
            U = triu(LU);
            testCase.verifyEqual(L*U, A, 'RelTol', 1e-9);
            testCase.verifyEqual(info, 0);
            % we do not the additional computations to enforce this!
            % testCase.verifyEqual(w/U, outw, 'RelTol', 1e-9);
            solver = str2func("mmex."+solvers{k} + "_solve");
            x = solver(A0, v, w, b, blocksize);
            testCase.verifyEqual(x, A\b, 'RelTol', 1e-9);
            solver = str2func("mmex."+solvers{k} + "_lsolve");
            x = solver(A0, v, w, b', blocksize);
            testCase.verifyEqual(x, b'/A, 'RelTol', 1e-9);
        end    
    end
end

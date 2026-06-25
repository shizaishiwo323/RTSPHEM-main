function tests = test_StepTransaction
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.rtmDir = rtmDir;
addpath(rtmDir);
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
end

function testRollbackReturnsOriginalSnapshots(testCase)
state = struct('time_s', 0, 'component_moles', [1 2]);
geometry = struct('water_volume_cm3', 1);
solverState = struct('dt_s', 10);
tx = rtm.driver.StepTransaction(state, geometry, solverState);

[rolledState, rolledGeometry, rolledSolverState] = tx.rollback();

verifyEqual(testCase, rolledState, state);
verifyEqual(testCase, rolledGeometry, geometry);
verifyEqual(testCase, rolledSolverState, solverState);
verifyEqual(testCase, tx.status(), 'rolled_back');
end

function testCommitReturnsCommittedSnapshots(testCase)
state = struct('time_s', 0);
geometry = struct('water_volume_cm3', 1);
solverState = struct();
tx = rtm.driver.StepTransaction(state, geometry, solverState);
newState = struct('time_s', 1);
newGeometry = struct('water_volume_cm3', 2);
newSolverState = struct('accepted_steps', 1);

[committedState, committedGeometry, committedSolverState] = ...
    tx.commit(newState, newGeometry, newSolverState);

verifyEqual(testCase, committedState, newState);
verifyEqual(testCase, committedGeometry, newGeometry);
verifyEqual(testCase, committedSolverState, newSolverState);
verifyEqual(testCase, tx.status(), 'committed');
end

function testRollbackAfterRollbackIsRejected(testCase)
tx = rtm.driver.StepTransaction(struct(), struct(), struct());
tx.rollback();

verifyError(testCase, @() tx.rollback(), 'RTSPHEM:Driver:TransactionClosed');
end

function testCommitAfterRollbackIsRejected(testCase)
tx = rtm.driver.StepTransaction(struct(), struct(), struct());
tx.rollback();

verifyError(testCase, @() tx.commit(struct(), struct(), struct()), ...
    'RTSPHEM:Driver:TransactionClosed');
end

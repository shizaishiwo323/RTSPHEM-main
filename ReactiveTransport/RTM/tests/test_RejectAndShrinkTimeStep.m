function tests = test_RejectAndShrinkTimeStep
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

function testRejectedStepShrinksDtAndRecordsReason(testCase)
solverState = struct('dt_s', 10, 'rejected_steps', 0, 'rejection_log', []);
options = struct('shrinkFactor', 0.5, 'minDt_s', 1, 'maxRetries', 3);

updated = rtm.driver.RejectAndShrinkTimeStep(solverState, "mass residual", options);

verifyEqual(testCase, updated.dt_s, 5);
verifyEqual(testCase, updated.rejected_steps, 1);
verifyEqual(testCase, updated.rejection_log(1).reason, "mass residual");
verifyEqual(testCase, updated.rejection_log(1).old_dt_s, 10);
verifyEqual(testCase, updated.rejection_log(1).new_dt_s, 5);
verifyFalse(testCase, updated.abort);
end

function testRetryLimitSetsAbortFlag(testCase)
solverState = struct('dt_s', 10, 'rejected_steps', 3, 'rejection_log', []);
options = struct('shrinkFactor', 0.5, 'minDt_s', 1, 'maxRetries', 3);

updated = rtm.driver.RejectAndShrinkTimeStep(solverState, "PHREEQC failed", options);

verifyTrue(testCase, updated.abort);
verifyEqual(testCase, updated.abort_reason, "maximum retries exceeded");
verifyEqual(testCase, updated.dt_s, 5);
end

function testMinimumDtSetsAbortFlag(testCase)
solverState = struct('dt_s', 1.5, 'rejected_steps', 0, 'rejection_log', []);
options = struct('shrinkFactor', 0.5, 'minDt_s', 1, 'maxRetries', 12);

updated = rtm.driver.RejectAndShrinkTimeStep(solverState, "geometry displacement", options);

verifyTrue(testCase, updated.abort);
verifyEqual(testCase, updated.abort_reason, "minimum dt reached");
verifyEqual(testCase, updated.dt_s, 0.75);
end

function testDefaultsMatchPlanFailureSettings(testCase)
solverState = struct('dt_s', 8);

updated = rtm.driver.RejectAndShrinkTimeStep(solverState, "default options");

verifyEqual(testCase, updated.dt_s, 4);
verifyEqual(testCase, updated.rejected_steps, 1);
verifyFalse(testCase, updated.abort);
end

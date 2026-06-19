function tests = test_StepperAppendTimeStepSize
tests = functiontests(localfunctions);
end

function testAppendTimeStepSizeExtendsStepper(testCase)
stepper = Stepper([0, 0.1, 0.3]);

verifyEqual(testCase, stepper.numsteps, 2);
verifyEqual(testCase, stepper.endtime, 0.3, 'AbsTol', 1e-12);

stepper.appendTimeStepSize(0.2);

verifyEqual(testCase, stepper.numsteps, 3);
verifyEqual(testCase, stepper.timepts(:)', [0, 0.1, 0.3, 0.5], 'AbsTol', 1e-12);
verifyEqual(testCase, stepper.endtime, 0.5, 'AbsTol', 1e-12);
end

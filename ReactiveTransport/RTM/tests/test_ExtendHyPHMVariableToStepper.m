function tests = test_ExtendHyPHMVariableToStepper
tests = functiontests(localfunctions);
end

function setupOnce(~)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rtmDir);
addpath(fullfile(fileparts(rtmDir), 'HyPHM', 'classes'));
end

function testCopyPreviousExtendsContentAndGrids(testCase)
stepper = Stepper([0, 1]);
fakeVariable = struct();
fakeVariable.stepper = stepper;
fakeVariable.content = {'old0'; 'old1'};
fakeVariable.grids = {'grid0'};

stepper.appendTimeStepSize(1);
fakeVariable = ExtendHyPHMVariableToStepper(fakeVariable, 'copy_previous');

verifyEqual(testCase, numel(fakeVariable.content), 3);
verifyEqual(testCase, fakeVariable.content{3}, 'old1');
verifyEqual(testCase, numel(fakeVariable.grids), 2);
verifyEqual(testCase, fakeVariable.grids{2}, 'grid0');
end

function testEmptyModeExtendsWithoutCopyingData(testCase)
stepper = Stepper([0, 1]);
fakeVariable = struct();
fakeVariable.stepper = stepper;
fakeVariable.content = {'old0'; 'old1'};
fakeVariable.grids = {'grid0'};

stepper.appendTimeStepSize(1);
fakeVariable = ExtendHyPHMVariableToStepper(fakeVariable, 'empty');

verifyEqual(testCase, numel(fakeVariable.content), 3);
verifyEmpty(testCase, fakeVariable.content{3});
verifyEqual(testCase, numel(fakeVariable.grids), 2);
end

function summary = RunQuasiSteadyGeometry(driver, totalTimeSeconds, options)
%RUNQUASISTEADYGEOMETRY Run geometry macro steps with rollback/retry.
%
% The helper provides a function-style entry point for benchmark orchestration
% while keeping retry state in the conservative driver.

if nargin < 3 || isempty(options)
    options = struct();
end
validateDriverMethod(driver, 'runQuasiSteadyGeometry');
summary = driver.runQuasiSteadyGeometry(totalTimeSeconds, options);
end

function validateDriverMethod(driver, methodName)
if isempty(driver) || ~ismethod(driver, methodName)
    error('RTSPHEM:Driver:InvalidDriver', ...
        'driver must provide %s.', methodName);
end
end

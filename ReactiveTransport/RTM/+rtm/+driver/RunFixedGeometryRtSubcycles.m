function summary = RunFixedGeometryRtSubcycles(driver, totalTimeSeconds, options)
%RUNFIXEDGEOMETRYRTSUBCYCLES Run RT subcycles while holding geometry fixed.
%
% This package-level helper is the functional entry point used by benchmark
% orchestration code. The implementation delegates to the conservative driver
% so there is a single source of time-step and retry behavior.

if nargin < 3 || isempty(options)
    options = struct();
end
validateDriverMethod(driver, 'runFixedGeometryRtSubcycles');
summary = driver.runFixedGeometryRtSubcycles(totalTimeSeconds, options);
end

function validateDriverMethod(driver, methodName)
if isempty(driver) || ~ismethod(driver, methodName)
    error('RTSPHEM:Driver:InvalidDriver', ...
        'driver must provide %s.', methodName);
end
end

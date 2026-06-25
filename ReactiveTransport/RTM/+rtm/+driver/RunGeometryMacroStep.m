function summary = RunGeometryMacroStep(driver, dtGeometrySeconds, options)
%RUNGEOMETRYMACROSTEP Run fixed-geometry RT, then move geometry once.
%
% This package-level helper keeps quasi-steady geometry orchestration out of
% benchmark scripts while delegating the actual conservative updates to the
% driver implementation.

if nargin < 3 || isempty(options)
    options = struct();
end
validateDriverMethod(driver, 'runGeometryMacroStep');
summary = driver.runGeometryMacroStep(dtGeometrySeconds, options);
end

function validateDriverMethod(driver, methodName)
if isempty(driver) || ~ismethod(driver, methodName)
    error('RTSPHEM:Driver:InvalidDriver', ...
        'driver must provide %s.', methodName);
end
end

function result = run_zhang_yoon_caco3_precipitation_benchmark(overrides)
% run_zhang_yoon_caco3_precipitation_benchmark
% Public Zhang 2010 / Yoon 2012 CaCO3 precipitation benchmark entrypoint.
%
% This entrypoint runs the Yoon Vaterite/Vm fixed-geometry micro-continuum
% path. The legacy signed Calcite level-set diagnostic is available through
% run_legacy_signed_calcite_levelset_diagnostic.

if nargin < 1
    overrides = struct();
end

moduleRoot = fileparts(mfilename('fullpath'));
rtmRoot = fileparts(moduleRoot);
reactiveRoot = fileparts(rtmRoot);

addpath(moduleRoot);
addpath(rtmRoot);
addpath(fullfile(rtmRoot, 'couplePhreeqc'));
addpath(genpath(fullfile(reactiveRoot, 'src')));
addpath(fullfile(reactiveRoot, 'HyPHM'));
addpath(fullfile(reactiveRoot, 'HyPHM', 'tools'));
addpath(genpath(fullfile(reactiveRoot, 'HyPHM', 'classes')));
addpath(genpath(fullfile(reactiveRoot, 'HyPHM', 'domains')));
addpath(genpath(fullfile(reactiveRoot, 'HyPHM', 'opt')));
addpath(genpath(fullfile(reactiveRoot, 'HyPHM', 'symbolic')));
addpath(moduleRoot, '-begin');
addpath(genpath(moduleRoot), '-begin');

specOverrides = getFieldOrDefault(overrides, 'spec', struct());
runOptions = getFieldOrDefault(overrides, 'options', struct());
spec = precip_ZhangYoonBenchmarkSpec(specOverrides);
result = precip_RunYoonFixedGeometryBenchmark(spec, runOptions);
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end

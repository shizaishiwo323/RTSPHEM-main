function result = run_legacy_signed_calcite_levelset_diagnostic(overrides)
% run_legacy_signed_calcite_levelset_diagnostic
% Legacy signed Calcite level-set diagnostic retained for comparison only.
%
% This path uses precip_ConfigureZhangYoonBenchmark and precip_PNM_beauty3.
% It is not the authoritative Zhang/Yoon Vaterite/Vm benchmark entrypoint.

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

cfg = precip_ConfigureZhangYoonBenchmark(overrides);
result = precip_PNM_beauty3(cfg);
end

function result = run_zhang_yoon_caco3_precipitation_benchmark(overrides)
% run_zhang_yoon_caco3_precipitation_benchmark
% Entry skeleton for the Zhang 2010 / Yoon 2012 CaCO3 precipitation benchmark.

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

cfg = precip_ConfigureZhangYoonBenchmark(overrides);
result = precip_PNM_beauty3(cfg);
end

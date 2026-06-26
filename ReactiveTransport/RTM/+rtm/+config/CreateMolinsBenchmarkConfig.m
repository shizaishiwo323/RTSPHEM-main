function config = CreateMolinsBenchmarkConfig(kind)
%CREATEMOLINSBENCHMARKCONFIG Build validated Molins benchmark configurations.

kind = lower(strrep(strtrim(char(kind)), '-', '_'));
config = baseConfig();

switch kind
    case {'integration_phreeqc', 'molins_geometry_phreeqc'}
        config.benchmark.name = 'molins_geometry_phreeqc_integration';
        config.benchmark.part = 'II';
        config.benchmark.strictBenchmark = false;
        config.benchmark.initializeFromPartI = true;
        config.chemistry.mode = 'external_tst_phreeqc';
        config.phreeqc.engine = 'iphreeqc_com';
        config.phreeqc.databasePolicy = 'exact_local';
        config.phreeqc.databaseName = 'phreeqc_rates.dat';
    otherwise
        error('RTSPHEM:Config:UnknownMolinsBenchmark', ...
            'Unknown Molins benchmark configuration: %s.', kind);
end

config = rtm.config.ValidateReactiveTransportConfig(config);
end

function config = baseConfig()
config = struct();
config.solverArchitecture = 'conservative_v2';
config.strictMassConservation = true;
config.chemistry = struct('mode', 'external_tst_phreeqc', ...
    'chargeAbsoluteTolerance_eq', 1e-8);
config.transport = struct('backend', 'cut_cell_fv');
config.time = struct();
config.time.mode = 'quasi_steady_geometry';
config.time.rt = struct('initialDt_s', 0.01, 'maxDt_s', 0.1, ...
    'advectiveCfl', 0.5, 'diffusiveNumber', 0.5, ...
    'maxReactantFraction', 0.1);
config.time.geometry = struct('boundaryCfl', 0.25, ...
    'maxMineralFraction', 0.2, 'maxDt_s', 60);
config.geometry = struct('molarVolume_cm3_mol', 1, ...
    'maxDisplacementOverH', 0.25, ...
    'solidAbsoluteTolerance_cm3', 1e-14, ...
    'solidRelativeTolerance', 1e-6);
config.cutCell = struct('minFluidFraction', 0.1, ...
    'smallCellMethod', 'disjoint_agglomeration');
config.mass = struct('absoluteTolerance_mol', 1e-14, ...
    'relativeTolerance', 1e-8, 'globalRelativeTolerance', 1e-6);
config.failure = struct('maxRetries', 12, 'shrinkFactor', 0.5, ...
    'minDt_s', 1e-8);
config.phreeqc = struct('engine', 'iphreeqc_com', 'databasePolicy', 'exact_local', ...
    'databaseName', 'phreeqc_rates.dat', 'persistSession', true, ...
    'useRunString', true);
config.benchmark = struct('enabled', true, 'name', '', 'part', '', ...
    'strictBenchmark', true, 'initializeFromPartI', false);
config.units = rtm.units.NormalizeRtmUnits();
end

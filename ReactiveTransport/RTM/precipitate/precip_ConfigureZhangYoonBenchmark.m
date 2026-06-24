function cfg = precip_ConfigureZhangYoonBenchmark(overrides)
% precip_ConfigureZhangYoonBenchmark - Configure the Zhang/Yoon CaCO3 case.

if nargin < 1
    overrides = struct();
end

moduleRoot = fileparts(mfilename('fullpath'));
rtmRoot = fileparts(moduleRoot);
reactiveRoot = fileparts(rtmRoot);
projectRoot = fileparts(reactiveRoot);

cfg = struct();
cfg.moduleRoot = moduleRoot;
cfg.rtmRoot = rtmRoot;
cfg.reactiveRoot = reactiveRoot;
cfg.projectRoot = projectRoot;
cfg.couplePhreeqcRoot = fullfile(rtmRoot, 'couplePhreeqc');

cfg.benchmarkName = 'zhang_yoon_caco3_25mM';
cfg.runName = [cfg.benchmarkName, '_', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))];
cfg.outputRoot = fullfile(moduleRoot, 'outputs', 'rtm_runs');
cfg.resultsDir = fullfile(cfg.outputRoot, cfg.runName);

cfg.reactionModel = 'phreeqc';
cfg.phreeqcRunGroup = 'phreeqc_database_calcite';
cfg.phreeqcRateLaw = 'phreeqc_database_calcite';
cfg.mineralEvolutionMode = 'signed_calcite_surface';

cfg.layoutType = 'zhang2010_micromodel_local';
cfg.lengthXAxis = 0.30;      % cm, local crop covering the inlet/first pores.
cfg.lengthYAxis = 0.30;      % cm, local crop transverse mixing region.
cfg.thickness = 0.002;       % cm, approximate 17-21 um Zhang/Yoon micromodel depth.
cfg.thicknessCm = cfg.thickness; % precip_PNM_beauty3 legacy interface field.
cfg.postDiameter = 0.03;     % cm, Zhang 2010 cylindrical post diameter.
cfg.poreBody = 0.018;        % cm, Zhang 2010 pore body width.
cfg.poreThroat = 0.004;      % cm, Zhang 2010 pore throat width.
cfg.targetPorosity = 0.39;   % -, Zhang 2010 micromodel porosity.
cfg.targetLengthYAxis = cfg.lengthYAxis;
cfg.targetAspectRatio = cfg.lengthXAxis / cfg.lengthYAxis;

cfg.flowBoundaryMode = 'split_left_inlet';
cfg.flowDirection = 'left_to_right';
cfg.splitInletY = 0.5 * cfg.lengthYAxis;
cfg.inletVelocity = 0.020833333;        % cm/s, Zhang 2010/Yoon 2012: 1.25 cm/min.
cfg.uD_cm_s = cfg.inletVelocity;        % cm/s, Zhang 2010/Yoon 2012 Darcy velocity.
cfg.diffusionCoefficient = 0.9e-5;      % cm2/s, Yoon 2012 aqueous diffusion value.
cfg.D_cm2_s = cfg.diffusionCoefficient; % cm2/s, Yoon 2012 notation.
cfg.molarVolume = 36.9;                 % cm3/mol, calcite effective molar volume.

cfg.benchmarkComparisonTimes_s = [13, 18, 118] * 60; % s, Zhang 2010/Yoon 2012 image times.
cfg.benchmarkComparisonTimesSource = 'Zhang2010_Yoon2012_image_comparison_times';
cfg.benchmarkFirstPoreXMaxCm = 0.078;       % cm, approximate local first-pore window.
cfg.benchmarkFirstThreePoresXMaxCm = 0.174; % cm, approximate first-three-pore window.
cfg.geometryFraming = 'approximate_local_cylindrical_post_layout_not_digitized_geometry';
cfg.endTime = 4 * 3600;
cfg.phreeqcExportEvery = 1;
cfg.exportEvery = 1;

cfg.inletA = struct();
cfg.inletA.name = 'CaCl2';
cfg.inletA.yRange = [0, cfg.splitInletY];
cfg.inletA.H_total = 10^(-6.1) * 1e-3;  % mol/cm3, Zhang 2010 CaCl2 inlet pH 6.1.
cfg.inletA.Ca_total = 25e-6;            % mol/cm3, Zhang 2010/Yoon 2012 25 mM CaCl2.
cfg.inletA.C_total = 0;                 % mol/cm3, Zhang 2010 CaCl2 inlet has no carbonate.
cfg.inletA.Na_total = 0;                % mol/cm3, Zhang 2010 CaCl2 inlet has no sodium.
cfg.inletA.Cl_total = 50e-6;            % mol/cm3, Zhang 2010/Yoon 2012 50 mM chloride.

cfg.inletB = struct();
cfg.inletB.name = 'Na2CO3';
cfg.inletB.yRange = [cfg.splitInletY, cfg.lengthYAxis];
cfg.inletB.H_total = 10^(-10.9) * 1e-3; % mol/cm3, Zhang 2010 Na2CO3 inlet pH 10.9.
cfg.inletB.Ca_total = 0;                % mol/cm3, Zhang 2010 Na2CO3 inlet has no calcium.
cfg.inletB.C_total = 25e-6;             % mol/cm3, Zhang 2010/Yoon 2012 25 mM carbonate.
cfg.inletB.Na_total = 50e-6;            % mol/cm3, Zhang 2010/Yoon 2012 50 mM sodium.
cfg.inletB.Cl_total = 0;                % mol/cm3, Zhang 2010 Na2CO3 inlet has no chloride.

% Keep scalar legacy inlets at the premixed stoichiometric average for
% uniform-inlet smoke tests and old scalar code paths.
cfg.initialHydrogenConcentration = 10^(-8.3) * 1e-3;
cfg.initialCalciumConcentration = 0;
cfg.initialCarbonConcentration = 0;
cfg.initialSodiumConcentration = 0;
cfg.initialChlorideConcentration = 0;
cfg.inletCalciumConcentration = 0.5 * cfg.inletA.Ca_total;
cfg.inletCarbonConcentration = 0.5 * cfg.inletB.C_total;
cfg.inletSodiumConcentration = 0.5 * cfg.inletB.Na_total;
cfg.inletChlorideConcentration = 0.5 * cfg.inletA.Cl_total;

cfg.precipitation = struct();
cfg.precipitation.enableSurfaceGrowth = true;
cfg.precipitation.enableHeterogeneousNucleation = false;
cfg.precipitation.enableHomogeneousNucleation = false;
cfg.precipitation.usePhreeqcSignedDelta = true;
cfg.precipitation.allowRedissolution = true;

cfg = mergeOverrides(cfg, overrides);
end

function cfg = mergeOverrides(cfg, overrides)
if ~isstruct(overrides)
    error('RTSPHEM:Precipitate:InvalidOverrides', 'Overrides must be a struct.');
end
hasInletAYRangeOverride = isNestedField(overrides, 'inletA', 'yRange');
hasInletBYRangeOverride = isNestedField(overrides, 'inletB', 'yRange');
fields = fieldnames(overrides);
for iField = 1:numel(fields)
    fieldName = fields{iField};
    if isfield(cfg, fieldName) && isstruct(cfg.(fieldName)) && isstruct(overrides.(fieldName))
        cfg.(fieldName) = mergeStructFields(cfg.(fieldName), overrides.(fieldName));
    else
        cfg.(fieldName) = overrides.(fieldName);
    end
end
if isfield(cfg, 'lengthYAxis') && (~isfield(overrides, 'splitInletY') || isempty(overrides.splitInletY))
    cfg.splitInletY = 0.5 * cfg.lengthYAxis;
end
if isfield(cfg, 'lengthYAxis') && isfield(cfg, 'splitInletY') && isfield(cfg, 'inletA') && isfield(cfg, 'inletB')
    if ~hasInletAYRangeOverride
    cfg.inletA.yRange = [0, cfg.splitInletY];
    end
    if ~hasInletBYRangeOverride
    cfg.inletB.yRange = [cfg.splitInletY, cfg.lengthYAxis];
    end
end
if isfield(cfg, 'thickness') && (~isfield(overrides, 'thicknessCm') || isempty(overrides.thicknessCm))
    cfg.thicknessCm = cfg.thickness;
end
if isfield(cfg, 'lengthXAxis') && isfield(cfg, 'lengthYAxis') && cfg.lengthYAxis > 0
    cfg.targetLengthYAxis = cfg.lengthYAxis;
    cfg.targetAspectRatio = cfg.lengthXAxis / cfg.lengthYAxis;
end
if isfield(cfg, 'inletA') && isfield(cfg, 'inletB')
    scalarFields = {'inletCalciumConcentration', 'inletCarbonConcentration', ...
        'inletSodiumConcentration', 'inletChlorideConcentration'};
    if ~any(isfield(overrides, scalarFields))
        cfg.inletCalciumConcentration = 0.5 * cfg.inletA.Ca_total;
        cfg.inletCarbonConcentration = 0.5 * cfg.inletB.C_total;
        cfg.inletSodiumConcentration = 0.5 * cfg.inletB.Na_total;
        cfg.inletChlorideConcentration = 0.5 * cfg.inletA.Cl_total;
    end
end
end

function merged = mergeStructFields(baseStruct, overrideStruct)
merged = baseStruct;
fields = fieldnames(overrideStruct);
for iField = 1:numel(fields)
    fieldName = fields{iField};
    merged.(fieldName) = overrideStruct.(fieldName);
end
end

function tf = isNestedField(parentStruct, fieldName, nestedFieldName)
tf = isstruct(parentStruct) && isfield(parentStruct, fieldName) && ...
    isstruct(parentStruct.(fieldName)) && isfield(parentStruct.(fieldName), nestedFieldName);
end

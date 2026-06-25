function spec = precip_ZhangYoonBenchmarkSpec(overrides)
% precip_ZhangYoonBenchmarkSpec - Source-traceable Zhang/Yoon benchmark constants.
%
% Inputs:
%   overrides - optional struct of scalar or nested fields to override.
%
% Output:
%   spec - benchmark constants for the Yoon Vaterite/Vm micro-continuum path.

if nargin < 1
    overrides = struct();
end
if ~isstruct(overrides)
    error('RTSPHEM:Precipitate:InvalidBenchmarkSpecOverrides', ...
        'Overrides must be a struct.');
end

spec = struct();
spec.name = 'zhang_yoon_caco3_25mM_yoon_case1';
spec.source = 'Zhang2010_Yoon2012';
spec.modelFamily = 'yoon_seeded_microcontinuum';
spec.precipitationMode = 'yoon_seeded_microcontinuum';
spec.deferredPrecipitationModes = {'deng_homogeneous_nucleation', ...
    'surface_growth'};
spec.deferredPrecipitationModePolicy = ...
    'fail_closed_until_yoon_quantitative_gate';
spec.mineralPhase = 'Vaterite';
spec.reactionRateLaw = 'yoon_chou_vaterite_explicit';
spec.solidStateVariable = 'Vm';
spec.chemistryBackend = 'yoon_equilibrium';
spec.componentNames = {'Ca_total', 'C_total', 'Na_total', 'Cl_total', 'Alkalinity'};

spec.lengthXAxis_cm = 0.30;
spec.lengthYAxis_cm = 0.30;
spec.thickness_cm = 0.002;
spec.numX = 60;
spec.numY = 60;
spec.splitInletY_cm = 0.5 * spec.lengthYAxis_cm;
spec.postDiameter_cm = 0.030;
spec.poreBodyWidth_cm = 0.018;
spec.poreThroatWidth_cm = 0.004;
spec.postPitchX_cm = spec.postDiameter_cm + spec.poreBodyWidth_cm;
spec.postPitchY_cm = spec.postDiameter_cm + spec.poreThroatWidth_cm;

spec.darcyVelocity_cm_s = 1.25 / 60;
spec.diffusionCoefficient_cm2_s = 9.0e-6;
spec.diffusionExponent = 2;
spec.blockedVmThreshold = 0.6;
spec.areaVmThreshold = 0.05;
spec.maxVmChangePerStep = 0.02;
spec.dissolutionFactor = 1;

spec.vateriteKsp_mol2_m6 = 1.832e-2;
spec.vateriteKsp_mol2_L2 = 1.832e-8;
spec.vateriteKsp = spec.vateriteKsp_mol2_L2;
spec.vateriteMolarVolume_cm3_mol = 37.47;
spec.temperature_C = 25;
spec.yoonRate = struct();
spec.yoonRate.k1 = 8.9e-5;
spec.yoonRate.k2 = 5.01e-8;
spec.yoonRate.k3 = 6.6e-11;
spec.yoonRate.units = 'mol_cm-2_s-1';
spec.yoonRate.source = 'Yoon2012_Chou_literature_locked';
spec.yoonRate.sourceDoi = '10.1029/2011WR011192';
spec.yoonRate.sourceEquation = 'Yoon2012 Eq.7; Chou1989';
spec.yoonRate.mineralPhase = 'Vaterite';
spec.yoonRate.sourceValuesVerified = true;

spec.inletA = makeInlet('CaCl2', [0, spec.splitInletY_cm], 6.1, ...
    25e-6, 0, 0, 50e-6);
spec.inletB = makeInlet('Na2CO3', [spec.splitInletY_cm, spec.lengthYAxis_cm], ...
    10.9, 0, 25e-6, 50e-6, 0);
spec.initial = makeInlet('DI_water', [0, spec.lengthYAxis_cm], 5.7, 0, 0, 0, 0);

spec = mergeOverrides(spec, overrides);
spec.dx_cm = spec.lengthXAxis_cm / spec.numX;
spec.dy_cm = spec.lengthYAxis_cm / spec.numY;
spec.cellVolume_cm3 = spec.dx_cm * spec.dy_cm * spec.thickness_cm;
spec.splitInletY_cm = resolveSplitY(spec);
spec.inletA.yRange_cm = resolveYRange(spec.inletA, [0, spec.splitInletY_cm]);
spec.inletB.yRange_cm = resolveYRange(spec.inletB, [spec.splitInletY_cm, spec.lengthYAxis_cm]);
end

function inlet = makeInlet(name, yRangeCm, pH, ca, carbon, na, cl)
inlet = struct();
inlet.name = name;
inlet.yRange_cm = yRangeCm;
inlet.pH = pH;
inlet.Ca_total = ca;
inlet.C_total = carbon;
inlet.Na_total = na;
inlet.Cl_total = cl;
inlet.Alkalinity = carbonateAlkalinity(carbon, pH);
end

function alkMolCm3 = carbonateAlkalinity(carbonMolCm3, pH)
carbonMolL = carbonMolCm3 * 1000;
h = 10^(-pH);
kw = 1e-14;
k1 = 10^(-6.35);
k2 = 10^(-10.33);
denom = h.^2 + k1 .* h + k1 .* k2;
alpha1 = k1 .* h ./ denom;
alpha2 = k1 .* k2 ./ denom;
alkMolL = carbonMolL .* (alpha1 + 2 .* alpha2) + kw ./ h - h;
alkMolCm3 = alkMolL / 1000;
end

function spec = mergeOverrides(spec, overrides)
fields = fieldnames(overrides);
for iField = 1:numel(fields)
    fieldName = fields{iField};
    if isfield(spec, fieldName) && isstruct(spec.(fieldName)) && ...
            isstruct(overrides.(fieldName))
        spec.(fieldName) = mergeStruct(spec.(fieldName), overrides.(fieldName));
    else
        spec.(fieldName) = overrides.(fieldName);
    end
end
end

function merged = mergeStruct(baseStruct, overrideStruct)
merged = baseStruct;
fields = fieldnames(overrideStruct);
for iField = 1:numel(fields)
    merged.(fields{iField}) = overrideStruct.(fields{iField});
end
end

function splitY = resolveSplitY(spec)
if isfield(spec, 'splitInletY_cm') && ~isempty(spec.splitInletY_cm)
    splitY = spec.splitInletY_cm;
else
    splitY = 0.5 * spec.lengthYAxis_cm;
end
end

function yRange = resolveYRange(inlet, fallback)
if isfield(inlet, 'yRange_cm') && numel(inlet.yRange_cm) == 2
    yRange = inlet.yRange_cm;
else
    yRange = fallback;
end
end

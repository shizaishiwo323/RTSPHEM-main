function suiteOptions = CreateMolinsDriverCaseOptions(kind, options)
%CREATEMOLINSDRIVERCASEOPTIONS Build convergence-suite options for Molins cases.

if nargin < 2 || isempty(options)
    options = struct();
end
kind = lower(strrep(strtrim(char(kind)), '-', '_'));
switch kind
    case {'integration_phreeqc', 'molins_geometry_phreeqc'}
        configKind = 'integration_phreeqc';
        suiteName = "molins_external_tst_phreeqc_driver";
        chemistryMode = 'external_tst_phreeqc';
    otherwise
        error('RTSPHEM:Benchmark:UnknownMolinsDriverCase', ...
            'Unknown Molins driver case kind: %s.', char(kind));
end

acceptanceCases = buildAcceptanceCases(options);
if isempty(acceptanceCases)
    refinementScales = getVectorOption(options, 'refinementScales', [5.4; 0.54; 0.054]);
else
    options.acceptanceCases = acceptanceCases;
    refinementScales = (numel(acceptanceCases):-1:1).';
end
caseOptions = struct();
caseOptions.totalTime_s = getScalarOption(options, 'totalTime_s', 1);
caseOptions.lightweightStepHistory = logical(getOption(options, ...
    'lightweightStepHistory', false));
caseOptions.useGeometryMacroLoop = logical(getOption(options, ...
    'useGeometryMacroLoop', ~caseOptions.lightweightStepHistory));
caseOptions.geometryMacroOptions = getOption(options, ...
    'geometryMacroOptions', struct());
if isfield(options, 'acceptanceMatrix') && ~isempty(options.acceptanceMatrix)
    caseOptions.acceptanceMatrix = options.acceptanceMatrix;
end
caseOptions.configFactory = @(scale, runInfo) molinsConfigFactory( ...
    scale, runInfo, configKind, chemistryMode, options);
caseOptions.stateFactory = @(scale, runInfo) molinsStateFactory( ...
    scale, runInfo, chemistryMode, options);
caseOptions.geometryFactory = @(scale, runInfo) molinsGeometryFactory( ...
    scale, runInfo, options);
caseOptions.connectivityFactory = @(scale, runInfo) molinsConnectivityFactory( ...
    scale, runInfo, options);
if ~isempty(acceptanceCases)
    caseOptions.acceptanceCases = acceptanceCases;
end

suiteOptions = struct();
suiteOptions.suiteName = suiteName;
suiteOptions.refinementScales = refinementScales(:);
if isempty(acceptanceCases)
    suiteOptions.runNames = runNamesFromScales(refinementScales);
else
    suiteOptions.runNames = runNamesFromAcceptanceCases(acceptanceCases);
end
suiteOptions.observableName = string(getOption(options, ...
    'observableName', 'final_porosity'));
suiteOptions.errorTolerance = getOption(options, 'errorTolerance', Inf);
suiteOptions.minObservedOrder = getOption(options, 'minObservedOrder', -Inf);
suiteOptions.referenceTargetValue = getOption(options, ...
    'referenceTargetValue', NaN);
suiteOptions.referenceRelativeTolerance = getOption(options, ...
    'referenceRelativeTolerance', Inf);
suiteOptions.referenceTargetRunNamePattern = getOption(options, ...
    'referenceTargetRunNamePattern', "");
suiteOptions.observableMaximumTolerance = getOption(options, ...
    'observableMaximumTolerance', Inf);
suiteOptions.refinementDimension = getOption(options, ...
    'refinementDimension', 'time');
suiteOptions.writePartialCheckpoint = logical(getOption(options, ...
    'writePartialCheckpoint', false));
suiteOptions.caseOptions = caseOptions;
if ~isempty(acceptanceCases)
    suiteOptions.acceptanceCases = acceptanceCases;
end
suiteOptions.runFunction = @(scale, runInfo) rtm.benchmark.RunDriverBenchmarkCase( ...
    scale, runInfo, caseOptions);
suiteOptions.observableFunction = observableFunctionFromName(suiteOptions.observableName);
if isfield(options, 'acceptanceMatrix') && ~isempty(options.acceptanceMatrix)
    suiteOptions.acceptanceMatrix = options.acceptanceMatrix;
end
if isfield(options, 'outputDir') && ~isempty(options.outputDir)
    suiteOptions.outputDir = options.outputDir;
end
end

function cfg = molinsConfigFactory(refinementScale, runInfo, configKind, chemistryMode, options)
caseSpec = acceptanceCaseForRun(refinementScale, runInfo, options);
if ~isempty(caseSpec)
    refinementScale = caseSpec.time_step_s;
end
cfg = rtm.config.CreateMolinsBenchmarkConfig(configKind);
cfg.time.rt.initialDt_s = refinementScale;
cfg.time.rt.maxDt_s = refinementScale;
cfg.time.rt.requestedDt_s = refinementScale;
cfg.chemistry.rate_constant_cm_s = getScalarOption(options, ...
    'rate_constant_cm_s', 0.1);
mesh = meshSpecForScale(refinementScale, options, caseSpec);
cfg.time.rt.maxReactantFraction = getFractionOption(options, ...
    'maxReactantFraction', cfg.time.rt.maxReactantFraction);
cfg.time.geometry.maxMineralFraction = getFractionOption(options, ...
    'maxMineralFraction', cfg.time.geometry.maxMineralFraction);
cfg.geometry.molarVolume_cm3_mol = getScalarOption(options, ...
    'molarVolume_cm3_mol', 1);
cfg.geometry.maxDisplacementOverH = getScalarOption(options, ...
    'maxDisplacementOverH', 0.25);
if ~isempty(mesh)
    cfg.benchmark.mesh = mesh;
    cfg.transport.options = molinsTransportOptions(mesh, options);
    if strcmp(cfg.transport.options.time_integration, 'explicit_euler')
        rtDt = molinsTransportCflDt(mesh, options);
        cfg.time.rt.initialDt_s = min(cfg.time.rt.initialDt_s, rtDt);
        cfg.time.rt.maxDt_s = min(cfg.time.rt.maxDt_s, rtDt);
    end
end
if strcmp(chemistryMode, 'external_tst_phreeqc')
    cfg.chemistry.options = struct();
    if isfield(cfg.transport, 'options') && isfield(cfg.transport.options, ...
            'boundary_face_cells')
        numBoundaryFaces = numel(cfg.transport.options.boundary_face_cells);
        cfg.transport.options.boundary_concentration_mol_cm3 = repmat([
            getScalarOption(options, 'inlet_ca_concentration_mol_cm3', 0), ...
            getScalarOption(options, 'inlet_c_concentration_mol_cm3', 0), ...
            getScalarOption(options, 'inlet_na_concentration_mol_cm3', 0), ...
            getScalarOption(options, 'inlet_cl_concentration_mol_cm3', 0)], ...
            numBoundaryFaces, 1);
    end
    cfg.phreeqc.databasePath = exactLocalPhreeqcDatabasePath( ...
        cfg.phreeqc.databaseName);
    cfg.chemistry.options.databasePath = cfg.phreeqc.databasePath;
    cfg.chemistry.options.h_mol_cm3 = getScalarOption(options, 'h_mol_cm3', 1e-7);
    cfg.chemistry.options.h_activity_mol_cm3 = getScalarOption(options, ...
        'h_activity_mol_cm3', cfg.chemistry.options.h_mol_cm3);
    cfg.chemistry.options.minReactionWaterVolumeCm3 = getScalarOption( ...
        options, 'minReactionWaterVolumeCm3', 1000);
    if isfield(options, 'phreeqcRunBatchFunction') && ~isempty(options.phreeqcRunBatchFunction)
        cfg.chemistry.options.runBatchFunction = options.phreeqcRunBatchFunction;
        cfg.phreeqc.engine = 'mock';
        cfg.phreeqc.databasePolicy = 'not_used';
        cfg.benchmark.enabled = false;
    end
end
end

function databasePath = exactLocalPhreeqcDatabasePath(databaseName)
benchmarkDir = fileparts(mfilename('fullpath'));
rtmDir = fileparts(fileparts(benchmarkDir));
databasePath = fullfile(rtmDir, 'phreeqc', 'database', databaseName);
if exist(databasePath, 'file') ~= 2
    error('RTSPHEM:Benchmark:MissingExactLocalDatabase', ...
        'Missing exact-local PHREEQC database: %s.', databasePath);
end
end

function state = molinsStateFactory(refinementScale, runInfo, chemistryMode, options)
caseSpec = acceptanceCaseForRun(refinementScale, runInfo, options);
if ~isempty(caseSpec)
    refinementScale = caseSpec.time_step_s;
end
numCells = molinsCellCount(refinementScale, options, caseSpec);
state = struct();
state.component_names = {'Ca', 'C', 'Na', 'Cl'};
totals = [
    getScalarOption(options, 'ca_moles', 0), ...
    getScalarOption(options, 'c_moles', 0), ...
    getScalarOption(options, 'na_moles', 1e-8), ...
    getScalarOption(options, 'cl_moles', 1e-8)];
mesh = meshSpecForScale(refinementScale, options, caseSpec);
if isempty(mesh)
    state.component_moles = distributeComponentTotals(totals, numCells);
else
    geometry = molinsGridGeometry(mesh, options);
    state.component_moles = distributeComponentTotalsOverWater( ...
        totals, geometry.water_volume_cm3(:));
end
state.mineral_names = {'Calcite'};
state.mineral_moles = molinsMineralMoles(refinementScale, options, caseSpec, numCells);
state.temperature_C = repmat(getScalarOption(options, 'temperature_C', 25), numCells, 1);
state.pressure_atm = repmat(getScalarOption(options, 'pressure_atm', 1), numCells, 1);
state.time_s = 0;
end

function connectivity = molinsConnectivityFactory(refinementScale, runInfo, options)
caseSpec = acceptanceCaseForRun(refinementScale, runInfo, options);
if ~isempty(caseSpec)
    refinementScale = caseSpec.time_step_s;
end
mesh = meshSpecForScale(refinementScale, options, caseSpec);
if isempty(mesh)
    connectivity = struct();
else
    connectivity = molinsGridConnectivity(mesh);
end
end

function values = molinsMineralMoles(refinementScale, options, caseSpec, numCells)
mesh = meshSpecForScale(refinementScale, options, caseSpec);
if isempty(mesh)
    calciteMoles = getScalarOption(options, 'calcite_moles', 1);
    values = distributeScalarMoles(calciteMoles, numCells);
    return;
end

geometry = molinsGridGeometry(mesh, options);
solidVolume = geometry.solid_volume_cm3(:);
molarVolume = getScalarOption(options, 'molarVolume_cm3_mol', 1);
if isfield(options, 'calcite_moles') && ~isempty(options.calcite_moles)
    calciteMoles = getScalarOption(options, 'calcite_moles', 1);
else
    calciteMoles = sum(solidVolume) ./ molarVolume;
end

totalSolidVolume = sum(solidVolume);
if totalSolidVolume > 0
    values = calciteMoles .* solidVolume ./ totalSolidVolume;
else
    values = distributeScalarMoles(calciteMoles, numCells);
end
end

function geometry = molinsGeometryFactory(refinementScale, runInfo, options)
caseSpec = acceptanceCaseForRun(refinementScale, runInfo, options);
if ~isempty(caseSpec)
    refinementScale = caseSpec.time_step_s;
end
mesh = meshSpecForScale(refinementScale, options, caseSpec);
if ~isempty(mesh)
    geometry = molinsGridGeometry(mesh, options);
    return;
end
waterVolume = getScalarOption(options, 'water_volume_cm3', 1);
solidVolume = getScalarOption(options, 'solid_volume_cm3', 1);
geometry = struct();
geometry.water_volume_cm3 = waterVolume;
geometry.solid_volume_cm3 = solidVolume;
geometry.cell_volume_cm3 = waterVolume + solidVolume;
geometry.fluid_fraction = waterVolume ./ max(geometry.cell_volume_cm3, eps);
geometry.cell_centroid_cm = [0 0];
geometry.interface_centroid_cm = [0 0];
geometry.interface_area_cm2 = getScalarOption(options, 'interface_area_cm2', 1);
geometry.interface_h_cm = getScalarOption(options, 'interface_h_cm', 1);
geometry.interface_normal = [1 0];
geometry.active_fluid_cell = true;
end

function mesh = meshSpecForScale(refinementScale, options, caseSpec)
if nargin >= 3 && ~isempty(caseSpec)
    if ~(isfield(caseSpec, 'mesh') && isstruct(caseSpec.mesh) && ...
            isfield(caseSpec.mesh, 'nx'))
        mesh = [];
        return;
    end
    mesh = caseSpec.mesh;
    return;
end
mesh = [];
if ~logical(getOption(options, 'useAcceptanceGrid', false))
    return;
end
if ~isfield(options, 'acceptanceMatrix') || isempty(options.acceptanceMatrix)
    return;
end
matrix = options.acceptanceMatrix;
if ~isfield(matrix, 'time_steps_s') || ~isfield(matrix, 'grid_resolutions')
    return;
end
timeSteps = matrix.time_steps_s(:);
matchIndex = find(abs(timeSteps - refinementScale) <= ...
    max(1e-12, 1e-12 .* abs(refinementScale)), 1);
if isempty(matchIndex)
    matchIndex = min(numel(matrix.grid_resolutions), numel(timeSteps));
end
labels = string(matrix.grid_resolutions(:));
matchIndex = min(matchIndex, numel(labels));
label = labels(matchIndex);
tokens = regexp(char(label), '^(\d+)x(\d+)$', 'tokens', 'once');
if isempty(tokens)
    error('RTSPHEM:Benchmark:InvalidAcceptanceGridResolution', ...
        'Invalid grid resolution label: %s.', char(label));
end
mesh = struct();
mesh.nx = str2double(tokens{1});
mesh.ny = str2double(tokens{2});
mesh.resolution_label = label;
mesh.refinement_scale_s = refinementScale;
mesh.domain_length_cm = getScalarOption(options, 'domain_length_cm', 0.1);
mesh.domain_height_cm = getScalarOption(options, 'domain_height_cm', 0.05);
mesh.circle_radius_cm = getScalarOption(options, 'circle_radius_cm', 0.01);
mesh.circle_center_cm = [ ...
    getScalarOption(options, 'circle_center_x_cm', 0.05), ...
    getScalarOption(options, 'circle_center_y_cm', 0.025)];
end

function numCells = molinsCellCount(refinementScale, options, caseSpec)
if nargin < 3
    caseSpec = [];
end
mesh = meshSpecForScale(refinementScale, options, caseSpec);
if isempty(mesh)
    numCells = 1;
else
    numCells = mesh.nx .* mesh.ny;
end
end

function values = distributeScalarMoles(totalMoles, numCells)
values = repmat(totalMoles ./ numCells, numCells, 1);
end

function values = distributeComponentTotals(totals, numCells)
values = repmat(totals ./ numCells, numCells, 1);
end

function values = distributeComponentTotalsOverWater(totals, waterVolume)
waterVolume = max(waterVolume(:), 0);
numCells = numel(waterVolume);
values = zeros(numCells, numel(totals));
totalWater = sum(waterVolume);
if totalWater <= 0
    return;
end
for iComponent = 1:numel(totals)
    values(:, iComponent) = totals(iComponent) .* waterVolume ./ totalWater;
end
end

function geometry = molinsGridGeometry(mesh, options)
nx = mesh.nx;
ny = mesh.ny;
dx = mesh.domain_length_cm ./ nx;
dy = mesh.domain_height_cm ./ ny;
[ix, iy] = ndgrid(1:nx, 1:ny);
centroidX = (ix(:) - 0.5) .* dx;
centroidY = (iy(:) - 0.5) .* dy;

cellArea = dx .* dy;
xMin = (ix(:) - 1) .* dx;
xMax = ix(:) .* dx;
yMin = (iy(:) - 1) .* dy;
yMax = iy(:) .* dy;

[solidVolume, waterVolume] = circleRectangleVolumes( ...
    xMin, xMax, yMin, yMax, cellArea, mesh, options);
[interfaceArea, interfaceCentroid, normal] = circleArcInterfaceMeasures( ...
    nx, ny, dx, dy, mesh, options);
cutCells = interfaceArea > 0 & solidVolume > 0 & waterVolume > 0;

cellVolume = repmat(cellArea, nx .* ny, 1);
geometry = struct();
geometry.water_volume_cm3 = waterVolume;
geometry.solid_volume_cm3 = solidVolume;
geometry.cell_volume_cm3 = cellVolume;
geometry.fluid_fraction = waterVolume ./ max(cellVolume, eps);
geometry.cell_centroid_cm = [centroidX, centroidY];
geometry.interface_centroid_cm = interfaceCentroid;
geometry.interface_normal = normal;
geometry.interface_area_cm2 = interfaceArea;
geometry.interface_h_cm = repmat(min(dx, dy), nx .* ny, 1);
geometry.interface_length_scale_cm = localInterfaceLengthScale( ...
    waterVolume, cellVolume, interfaceArea);
geometry.active_fluid_cell = waterVolume > 0;
geometry.cut_cell = cutCells;
geometry.mesh_resolution = [nx ny];
geometry.mesh_resolution_label = mesh.resolution_label;
geometry.initial_surface_area_cm2 = sum(interfaceArea);
geometry.initial_solid_volume_cm3 = sum(solidVolume);
geometry.domain_length_cm = mesh.domain_length_cm;
geometry.domain_height_cm = mesh.domain_height_cm;
end

function lengthScale = localInterfaceLengthScale(waterVolume, cellVolume, interfaceArea)
lengthScale = zeros(size(interfaceArea));
activeInterface = interfaceArea > 0;
lengthScale(activeInterface) = min( ...
    sqrt(waterVolume(activeInterface) ./ max(interfaceArea(activeInterface), eps)), ...
    sqrt(cellVolume(activeInterface)));
end

function [solidVolume, waterVolume] = circleRectangleVolumes( ...
        xMin, xMax, yMin, yMax, cellArea, mesh, options)
numCells = numel(xMin);
solidVolume = zeros(numCells, 1);
center = mesh.circle_center_cm;
radius = mesh.circle_radius_cm;

closestDx = max(max(xMin - center(1), 0), center(1) - xMax);
closestDy = max(max(yMin - center(2), 0), center(2) - yMax);
closestDistance = hypot(closestDx, closestDy);

cornerDistance = [ ...
    hypot(xMin - center(1), yMin - center(2)), ...
    hypot(xMin - center(1), yMax - center(2)), ...
    hypot(xMax - center(1), yMin - center(2)), ...
    hypot(xMax - center(1), yMax - center(2))];
farthestDistance = max(cornerDistance, [], 2);

fullSolid = farthestDistance <= radius;
fullFluid = closestDistance >= radius;
cutCells = ~(fullSolid | fullFluid);

solidVolume(fullSolid) = cellArea;
sampleCount = max(1, round(getScalarOption(options, ...
    'circle_cut_subsamples', 128)));
offsets = ((1:sampleCount) - 0.5) ./ sampleCount;
[ox, oy] = ndgrid(offsets, offsets);
ox = ox(:);
oy = oy(:);
for iCell = find(cutCells(:)).'
    sampleX = xMin(iCell) + ox .* (xMax(iCell) - xMin(iCell));
    sampleY = yMin(iCell) + oy .* (yMax(iCell) - yMin(iCell));
    insideFraction = mean(hypot(sampleX - center(1), ...
        sampleY - center(2)) <= radius);
    solidVolume(iCell) = insideFraction .* cellArea;
end

solidVolume = min(max(solidVolume, 0), cellArea);
waterVolume = cellArea - solidVolume;
end

function [interfaceArea, interfaceCentroid, normal] = circleArcInterfaceMeasures( ...
        nx, ny, dx, dy, mesh, options)
numCells = nx .* ny;
interfaceArea = zeros(numCells, 1);
interfaceCentroid = nan(numCells, 2);
normal = nan(numCells, 2);

radius = mesh.circle_radius_cm;
center = mesh.circle_center_cm;
sampleCount = max(4096, round(getScalarOption(options, ...
    'interface_arc_samples', 32768)));
theta = (0:sampleCount - 1).' .* (2 .* pi ./ sampleCount);
arcX = center(1) + radius .* cos(theta);
arcY = center(2) + radius .* sin(theta);

ix = min(max(floor(arcX ./ dx) + 1, 1), nx);
iy = min(max(floor(arcY ./ dy) + 1, 1), ny);
cellId = ix + (iy - 1) .* nx;
arcWeight = 2 .* pi .* radius ./ sampleCount;

sumX = accumarray(cellId, arcX, [numCells, 1], @sum, 0);
sumY = accumarray(cellId, arcY, [numCells, 1], @sum, 0);
counts = accumarray(cellId, 1, [numCells, 1], @sum, 0);
interfaceArea = counts .* arcWeight;

active = counts > 0;
meanX = zeros(numCells, 1);
meanY = zeros(numCells, 1);
meanX(active) = sumX(active) ./ counts(active);
meanY(active) = sumY(active) ./ counts(active);
radial = [meanX - center(1), meanY - center(2)];
radialLength = hypot(radial(:, 1), radial(:, 2));
project = active & radialLength > 0;
interfaceCentroid(project, :) = center + ...
    radius .* radial(project, :) ./ radialLength(project);
normal(project, :) = radial(project, :) ./ radialLength(project);
axisFallback = active & ~project;
interfaceCentroid(axisFallback, :) = repmat(center + [radius 0], ...
    nnz(axisFallback), 1);
normal(axisFallback, :) = repmat([1 0], nnz(axisFallback), 1);

interfaceAreaTotal = getScalarOption(options, 'interface_area_total_cm2', ...
    2 .* pi .* radius);
if sum(interfaceArea) > 0
    interfaceArea = interfaceArea .* (interfaceAreaTotal ./ sum(interfaceArea));
end
end

function connectivity = molinsGridConnectivity(mesh)
nx = mesh.nx;
ny = mesh.ny;
numCells = nx .* ny;
neighbors = repmat({zeros(0, 1)}, numCells, 1);
for iy = 1:ny
    for ix = 1:nx
        cellId = gridCellId(ix, iy, nx);
        ids = zeros(0, 1);
        if ix > 1
            ids(end + 1, 1) = gridCellId(ix - 1, iy, nx); %#ok<AGROW>
        end
        if ix < nx
            ids(end + 1, 1) = gridCellId(ix + 1, iy, nx); %#ok<AGROW>
        end
        if iy > 1
            ids(end + 1, 1) = gridCellId(ix, iy - 1, nx); %#ok<AGROW>
        end
        if iy < ny
            ids(end + 1, 1) = gridCellId(ix, iy + 1, nx); %#ok<AGROW>
        end
        neighbors{cellId} = ids;
    end
end
connectivity = struct('cell_neighbors', {neighbors});
end

function transportOptions = molinsTransportOptions(mesh, options)
geometry = molinsGridGeometry(mesh, options);
nx = mesh.nx;
ny = mesh.ny;
dx = mesh.domain_length_cm ./ nx;
dy = mesh.domain_height_cm ./ ny;
active = geometry.water_volume_cm3(:) > 0;

faceCells = zeros(0, 2);
faceArea = zeros(0, 1);
faceDistance = zeros(0, 1);
faceVelocity = zeros(0, 1);
inletVelocity = getScalarOption(options, 'inlet_velocity_cm_s', 0.12);
verticalOpenTotals = molinsVerticalOpenTotals(mesh, active, nx, ny, dx, dy);

for iy = 1:ny
    for ix = 1:nx
        leftCell = gridCellId(ix, iy, nx);
        if ix < nx
            rightCell = gridCellId(ix + 1, iy, nx);
            if active(leftCell) && active(rightCell)
                yLower = (iy - 1) .* dy;
                yUpper = iy .* dy;
                xFace = ix .* dx;
                openArea = circleVerticalOpenLength(xFace, yLower, yUpper, mesh);
                if openArea > 0
                    columnVelocity = inletVelocity .* mesh.domain_height_cm ./ ...
                        max(verticalOpenTotals(ix), eps);
                    faceCells(end + 1, :) = [leftCell rightCell]; %#ok<AGROW>
                    faceArea(end + 1, 1) = openArea; %#ok<AGROW>
                    faceDistance(end + 1, 1) = dx; %#ok<AGROW>
                    faceVelocity(end + 1, 1) = columnVelocity; %#ok<AGROW>
                end
            end
        end
        if iy < ny
            topCell = gridCellId(ix, iy + 1, nx);
            if active(leftCell) && active(topCell)
                xLower = (ix - 1) .* dx;
                xUpper = ix .* dx;
                yFace = iy .* dy;
                openArea = circleHorizontalOpenLength(yFace, xLower, xUpper, mesh);
                if openArea > 0
                    verticalVelocity = molinsHorizontalFaceVelocity( ...
                        yFace, xLower, xUpper, mesh, inletVelocity, options);
                    faceCells(end + 1, :) = [leftCell topCell]; %#ok<AGROW>
                    faceArea(end + 1, 1) = openArea; %#ok<AGROW>
                    faceDistance(end + 1, 1) = dy; %#ok<AGROW>
                    faceVelocity(end + 1, 1) = verticalVelocity; %#ok<AGROW>
                end
            end
        end
    end
end

leftBoundaryCells = zeros(0, 1);
leftBoundaryArea = zeros(0, 1);
rightBoundaryCells = zeros(0, 1);
rightBoundaryArea = zeros(0, 1);
for iy = 1:ny
    cellId = gridCellId(1, iy, nx);
    if active(cellId)
        leftBoundaryCells(end + 1, 1) = cellId; %#ok<AGROW>
        yLower = (iy - 1) .* dy;
        yUpper = iy .* dy;
        leftBoundaryArea(end + 1, 1) = ...
            circleVerticalOpenLength(0, yLower, yUpper, mesh); %#ok<AGROW>
    end

    rightCellId = gridCellId(nx, iy, nx);
    if active(rightCellId)
        rightBoundaryCells(end + 1, 1) = rightCellId; %#ok<AGROW>
        yLower = (iy - 1) .* dy;
        yUpper = iy .* dy;
        rightBoundaryArea(end + 1, 1) = ...
            circleVerticalOpenLength(mesh.domain_length_cm, yLower, yUpper, mesh); %#ok<AGROW>
    end
end

boundaryCells = [leftBoundaryCells; rightBoundaryCells];
boundaryArea = [leftBoundaryArea; rightBoundaryArea];
boundaryType = [repmat("dirichlet", numel(leftBoundaryCells), 1); ...
    repmat("outflow", numel(rightBoundaryCells), 1)];

transportOptions = struct();
transportOptions.time_integration = char(getOption(options, ...
    'transport_time_integration', 'implicit_euler'));
transportOptions.internal_face_cells = faceCells;
transportOptions.internal_face_area_cm2 = faceArea;
transportOptions.internal_face_distance_cm = faceDistance;
transportOptions.internal_face_velocity_cm_s = faceVelocity;
transportOptions.boundary_face_cells = boundaryCells;
transportOptions.boundary_face_area_cm2 = boundaryArea;
transportOptions.boundary_face_distance_cm = repmat(dx ./ 2, numel(boundaryCells), 1);
transportOptions.boundary_face_velocity_cm_s = repmat(inletVelocity, ...
    numel(boundaryCells), 1);
transportOptions.boundary_type = boundaryType;
transportOptions.boundary_concentration_mol_cm3 = repmat( ...
    getScalarOption(options, 'inlet_h_concentration_mol_cm3', 1.255e-6), ...
    numel(boundaryCells), 1);
transportOptions.diffusion_coefficient_cm2_s = getScalarOption(options, ...
    'diffusion_coefficient_cm2_s', 1e-5);
end

function totals = molinsVerticalOpenTotals(mesh, active, nx, ny, dx, dy)
totals = zeros(nx - 1, 1);
for ix = 1:nx - 1
    xFace = ix .* dx;
    for iy = 1:ny
        leftCell = gridCellId(ix, iy, nx);
        rightCell = gridCellId(ix + 1, iy, nx);
        if ~(active(leftCell) && active(rightCell))
            continue;
        end
        yLower = (iy - 1) .* dy;
        yUpper = iy .* dy;
        totals(ix) = totals(ix) + ...
            circleVerticalOpenLength(xFace, yLower, yUpper, mesh);
    end
end
end

function openLength = circleVerticalOpenLength(xValue, yLower, yUpper, mesh)
radius = mesh.circle_radius_cm;
center = mesh.circle_center_cm;
solidHalfChordSquared = radius.^2 - (xValue - center(1)).^2;
solidLength = 0;
if solidHalfChordSquared > 0
    halfChord = sqrt(solidHalfChordSquared);
    solidLower = center(2) - halfChord;
    solidUpper = center(2) + halfChord;
    solidLength = max(0, min(yUpper, solidUpper) - max(yLower, solidLower));
end
openLength = max(0, (yUpper - yLower) - solidLength);
end

function openLength = circleHorizontalOpenLength(yValue, xLower, xUpper, mesh)
radius = mesh.circle_radius_cm;
center = mesh.circle_center_cm;
solidHalfChordSquared = radius.^2 - (yValue - center(2)).^2;
solidLength = 0;
if solidHalfChordSquared > 0
    halfChord = sqrt(solidHalfChordSquared);
    solidLower = center(1) - halfChord;
    solidUpper = center(1) + halfChord;
    solidLength = max(0, min(xUpper, solidUpper) - max(xLower, solidLower));
end
openLength = max(0, (xUpper - xLower) - solidLength);
end

function velocity = molinsHorizontalFaceVelocity( ...
        yValue, xLower, xUpper, mesh, inletVelocity, options)
if ~logical(getOption(options, 'usePotentialFlowVelocity', false))
    velocity = 0;
    return;
end
intervals = circleHorizontalOpenIntervals(yValue, xLower, xUpper, mesh);
if isempty(intervals)
    velocity = 0;
    return;
end
weightedVelocity = 0;
openLength = 0;
for iInterval = 1:size(intervals, 1)
    intervalLength = intervals(iInterval, 2) - intervals(iInterval, 1);
    if intervalLength <= 0
        continue;
    end
    xMid = mean(intervals(iInterval, :));
    [~, velocityY] = potentialCylinderVelocity(xMid, yValue, mesh, inletVelocity);
    weightedVelocity = weightedVelocity + velocityY .* intervalLength;
    openLength = openLength + intervalLength;
end
if openLength > 0
    velocity = weightedVelocity ./ openLength;
else
    velocity = 0;
end
end

function intervals = circleHorizontalOpenIntervals(yValue, xLower, xUpper, mesh)
radius = mesh.circle_radius_cm;
center = mesh.circle_center_cm;
solidHalfChordSquared = radius.^2 - (yValue - center(2)).^2;
if solidHalfChordSquared <= 0
    intervals = [xLower, xUpper];
    return;
end
halfChord = sqrt(solidHalfChordSquared);
solidLower = max(xLower, center(1) - halfChord);
solidUpper = min(xUpper, center(1) + halfChord);
if solidUpper <= xLower || solidLower >= xUpper || solidUpper <= solidLower
    intervals = [xLower, xUpper];
    return;
end
intervals = zeros(0, 2);
if solidLower > xLower
    intervals(end + 1, :) = [xLower, solidLower]; %#ok<AGROW>
end
if solidUpper < xUpper
    intervals(end + 1, :) = [solidUpper, xUpper]; %#ok<AGROW>
end
end

function [velocityX, velocityY] = potentialCylinderVelocity( ...
        xValue, yValue, mesh, inletVelocity)
relativeX = xValue - mesh.circle_center_cm(1);
relativeY = yValue - mesh.circle_center_cm(2);
radiusSquared = relativeX.^2 + relativeY.^2;
cylinderRadius = mesh.circle_radius_cm;
if radiusSquared <= cylinderRadius.^2
    velocityX = 0;
    velocityY = 0;
    return;
end
radiusFourth = radiusSquared.^2;
velocityX = inletVelocity .* ...
    (1 - cylinderRadius.^2 .* (relativeX.^2 - relativeY.^2) ./ radiusFourth);
velocityY = -2 .* inletVelocity .* cylinderRadius.^2 .* ...
    relativeX .* relativeY ./ radiusFourth;
end

function dtSeconds = molinsTransportCflDt(mesh, options)
dx = mesh.domain_length_cm ./ mesh.nx;
dy = mesh.domain_height_cm ./ mesh.ny;
h = min(dx, dy);
velocity = getScalarOption(options, 'inlet_velocity_cm_s', 0.12);
diffusion = getScalarOption(options, 'diffusion_coefficient_cm2_s', 1e-5);
advectiveCfl = getScalarOption(options, 'advective_cfl', 0.25);
diffusiveNumber = getScalarOption(options, 'diffusive_number', 0.25);
dtAdv = Inf;
if velocity > 0
    dtAdv = advectiveCfl .* h ./ velocity;
end
dtDiff = Inf;
if diffusion > 0
    dtDiff = diffusiveNumber .* h.^2 ./ diffusion;
end
dtSeconds = max(min(dtAdv, dtDiff), 1e-8);
end

function cellId = gridCellId(ix, iy, nx)
cellId = ix + (iy - 1) .* nx;
end

function handle = observableFunctionFromName(observableName)
switch char(observableName)
    case 'final_porosity'
        handle = @(summary) summary.final_porosity;
    case 'mineral_dissolved_moles'
        handle = @(summary) summary.mineral_dissolved_moles;
    case 'final_mineral_moles'
        handle = @(summary) summary.final_mineral_moles;
    otherwise
        handle = @(summary) summary.(char(observableName));
end
end

function names = runNamesFromScales(scales)
names = strings(numel(scales), 1);
for iScale = 1:numel(scales)
    text = regexprep(sprintf('%.15g', scales(iScale)), '[^0-9A-Za-z]', 'p');
    names(iScale) = "dt_" + text;
end
end

function names = runNamesFromAcceptanceCases(cases)
names = strings(numel(cases), 1);
for iCase = 1:numel(cases)
    names(iCase) = cases(iCase).name;
end
end

function cases = buildAcceptanceCases(options)
cases = struct([]);
if ~isfield(options, 'acceptanceMatrix') || isempty(options.acceptanceMatrix)
    return;
end
matrix = options.acceptanceMatrix;
if ~isfield(matrix, 'time_steps_s') || ~isfield(matrix, 'grid_resolutions')
    return;
end
timeSteps = matrix.time_steps_s(:);
gridLabels = string(matrix.grid_resolutions(:));
if any(~isfinite(timeSteps)) || any(timeSteps <= 0)
    error('RTSPHEM:Benchmark:InvalidAcceptanceMatrix', ...
        'acceptanceMatrix.time_steps_s must contain positive finite values.');
end
if isempty(gridLabels)
    error('RTSPHEM:Benchmark:InvalidAcceptanceMatrix', ...
        'acceptanceMatrix.grid_resolutions must not be empty.');
end
cases = repmat(emptyAcceptanceCase(), numel(timeSteps) * numel(gridLabels), 1);
iCase = 0;
for iGrid = 1:numel(gridLabels)
    if logical(getOption(options, 'useAcceptanceGrid', false))
        mesh = meshFromGridLabel(gridLabels(iGrid), options);
    else
        mesh = struct();
    end
    for iTime = 1:numel(timeSteps)
        iCase = iCase + 1;
        cases(iCase).time_step_s = timeSteps(iTime);
        cases(iCase).grid_resolution = gridLabels(iGrid);
        cases(iCase).grid_spacing_cm = gridSpacingFromLabel(gridLabels(iGrid), ...
            options);
        cases(iCase).mesh = mesh;
        cases(iCase).name = "grid_" + gridLabels(iGrid) + "_dt_" + ...
            sanitizeRunToken(timeSteps(iTime));
    end
end
end

function value = emptyAcceptanceCase()
value = struct();
value.time_step_s = NaN;
value.grid_resolution = "";
value.grid_spacing_cm = NaN;
value.mesh = struct();
value.name = "";
end

function value = gridSpacingFromLabel(label, options)
mesh = meshFromGridLabel(label, options);
value = min(mesh.domain_length_cm ./ mesh.nx, mesh.domain_height_cm ./ mesh.ny);
end

function caseSpec = acceptanceCaseForRun(refinementScale, runInfo, options)
caseSpec = [];
if ~isfield(options, 'acceptanceCases') || isempty(options.acceptanceCases)
    return;
end
cases = options.acceptanceCases;
index = NaN;
if isstruct(runInfo) && isfield(runInfo, 'index') && ~isempty(runInfo.index)
    index = runInfo.index;
elseif isscalar(refinementScale) && isfinite(refinementScale) && ...
        refinementScale == round(refinementScale)
    index = refinementScale;
end
if ~(isscalar(index) && isfinite(index) && index >= 1 && ...
        index <= numel(cases)) && isscalar(refinementScale) && ...
        isfinite(refinementScale)
    timeSteps = [cases.time_step_s].';
    tolerance = max(1e-12, 1e-12 .* max(abs(refinementScale), 1));
    index = find(abs(timeSteps - refinementScale) <= tolerance, 1, 'first');
end
if ~(isscalar(index) && isfinite(index) && index >= 1 && ...
        index <= numel(cases))
    error('RTSPHEM:Benchmark:InvalidAcceptanceCaseIndex', ...
        'Acceptance case index is out of range.');
end
caseSpec = cases(index);
end

function mesh = meshFromGridLabel(label, options)
tokens = regexp(char(label), '^(\d+)x(\d+)$', 'tokens', 'once');
if isempty(tokens)
    error('RTSPHEM:Benchmark:InvalidAcceptanceGridResolution', ...
        'Invalid grid resolution label: %s.', char(label));
end
mesh = struct();
mesh.nx = str2double(tokens{1});
mesh.ny = str2double(tokens{2});
mesh.resolution_label = string(label);
mesh.domain_length_cm = getScalarOption(options, 'domain_length_cm', 0.1);
mesh.domain_height_cm = getScalarOption(options, 'domain_height_cm', 0.05);
mesh.circle_radius_cm = getScalarOption(options, 'circle_radius_cm', 0.01);
mesh.circle_center_cm = [ ...
    getScalarOption(options, 'circle_center_x_cm', 0.05), ...
    getScalarOption(options, 'circle_center_y_cm', 0.025)];
end

function token = sanitizeRunToken(value)
token = regexprep(sprintf('%.15g', value), '[^0-9A-Za-z]', 'p');
end

function values = getVectorOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    values = options.(fieldName);
else
    values = defaultValue;
end
values = values(:);
if any(~isfinite(values)) || any(values <= 0)
    error('RTSPHEM:Benchmark:InvalidMolinsDriverCaseOption', ...
        'options.%s must contain positive finite values.', fieldName);
end
end

function value = getScalarOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
if ~(isscalar(value) && isfinite(value))
    error('RTSPHEM:Benchmark:InvalidMolinsDriverCaseOption', ...
        'options.%s must be a finite scalar.', fieldName);
end
end

function value = getFractionOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
if ~(isscalar(value) && (isfinite(value) || isequal(value, Inf)) && value >= 0)
    error('RTSPHEM:Benchmark:InvalidMolinsDriverCaseOption', ...
        'options.%s must be a nonnegative finite scalar or Inf.', fieldName);
end
end

function value = getIntegerOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
if ~(isscalar(value) && isfinite(value) && value >= 0 && value == round(value))
    error('RTSPHEM:Benchmark:InvalidMolinsDriverCaseOption', ...
        'options.%s must be a nonnegative integer scalar.', fieldName);
end
end

function value = getOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

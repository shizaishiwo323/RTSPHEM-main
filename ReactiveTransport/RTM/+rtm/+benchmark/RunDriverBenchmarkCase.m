function summary = RunDriverBenchmarkCase(refinementScale, runInfo, options)
%RUNDRIVERBENCHMARKCASE Execute one ReactiveTransportDriver benchmark case.

if nargin < 2 || isempty(runInfo)
    runInfo = struct();
end
if nargin < 3 || isempty(options)
    options = struct();
end
validateFactories(options);

cfg = options.configFactory(refinementScale, runInfo);
state = options.stateFactory(refinementScale, runInfo);
geometry = options.geometryFactory(refinementScale, runInfo);
connectivity = buildConnectivity(options, refinementScale, runInfo);
totalTimeSeconds = getScalarOption(options, 'totalTime_s', 1);

initialState = state;
initialGeometry = geometry;
try
    steadyInfo = struct();
    if shouldInitializeFromPartI(cfg)
        [state, geometry, steadyInfo] = initializeFromPartISteadyState( ...
            cfg, state, geometry, connectivity, options);
    end
    driver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, connectivity);
    if shouldUseLightweightStepHistory(options)
        driverSummary = runRtSubcyclesLightweight(driver, totalTimeSeconds, ...
            initialState, initialGeometry, cfg, options);
    else
        driverSummary = driver.runRtSubcycles(totalTimeSeconds);
    end
    summary = enrichSummary(driverSummary, initialState, initialGeometry, ...
        cfg, refinementScale, runInfo, options);
    if ~isempty(fieldnames(steadyInfo))
        summary.steady_initializer = steadyInfo;
    end
catch err
    summary = rejectedSummary(err, initialState, initialGeometry, cfg, ...
        refinementScale, runInfo, options);
end
end

function tf = shouldInitializeFromPartI(cfg)
tf = logical(getNestedField(cfg, {'benchmark', 'initializeFromPartI'}, false));
end

function [initializedState, initializedGeometry, steadyInfo] = ...
    initializeFromPartISteadyState(cfg, state, geometry, connectivity, options)
steadyDriver = rtm.driver.ReactiveTransportDriver(cfg, state, geometry, connectivity);
steadyOptions = getFieldOrDefault(options, 'steadyOptions', struct());
steadyInfo = rtm.driver.RunToReactiveSteadyState(steadyDriver, steadyOptions);

initializedState = steadyInfo.state;
initializedState = remapInitializerStateToGeometry( ...
    initializedState, steadyInfo.geometry, geometry, state.mineral_moles);
initializedState.time_s = state.time_s;
initializedGeometry = geometry;

if isfield(steadyInfo, 'geometry') && isstruct(steadyInfo.geometry)
    steadyInfo.initializer_remap = initializerRemapSummary( ...
        steadyInfo.state, steadyInfo.geometry, initializedState, geometry);
end
steadyInfo.state = stripHeavyField(steadyInfo.state);
steadyInfo.geometry = stripHeavyField(steadyInfo.geometry);
steadyInfo.summaries = struct([]);
end

function initializedState = remapInitializerStateToGeometry( ...
        steadyState, steadyGeometry, targetGeometry, targetMineralMoles)
options = struct('new_mineral_moles', targetMineralMoles);
[initializedState, ledger] = rtm.geometry.ConservativeRemap( ...
    steadyState, steadyGeometry, targetGeometry, options);
unassigned = ledger.unassigned_component_moles;
if any(abs(unassigned) > 0)
    initializedState.component_moles = redistributeUnassignedComponents( ...
        initializedState.component_moles, targetGeometry, unassigned);
end
initializedState.mineral_moles = targetMineralMoles;
end

function componentMoles = redistributeUnassignedComponents( ...
        componentMoles, geometry, unassignedComponentMoles)
waterVolume = geometry.water_volume_cm3(:);
active = waterVolume > 0;
if ~any(active)
    if any(abs(unassignedComponentMoles) > 1e-18)
        error('RTSPHEM:Benchmark:InitializerRemapNoTargetWater', ...
            'Initializer remap has unassigned component moles but no target water cells.');
    end
    return;
end
weights = waterVolume(active) ./ sum(waterVolume(active));
for iComponent = 1:numel(unassignedComponentMoles)
    componentMoles(active, iComponent) = componentMoles(active, iComponent) + ...
        weights .* unassignedComponentMoles(iComponent);
end
componentMoles(~isfinite(componentMoles)) = 0;
end

function summary = initializerRemapSummary( ...
        steadyState, steadyGeometry, initializedState, targetGeometry)
summary = struct();
summary.initial_component_moles_total = sum(steadyState.component_moles, 1);
summary.final_component_moles_total = sum(initializedState.component_moles, 1);
summary.component_residual_moles = ...
    summary.final_component_moles_total - summary.initial_component_moles_total;
targetDry = targetGeometry.water_volume_cm3(:) <= 0;
summary.dry_cell_component_moles_total = ...
    sum(initializedState.component_moles(targetDry, :), 1);
summary.steady_water_volume_cm3 = sum(steadyGeometry.water_volume_cm3(:), 'omitnan');
summary.target_water_volume_cm3 = sum(targetGeometry.water_volume_cm3(:), 'omitnan');
end

function value = stripHeavyField(value)
% Keep the initializer manifest lightweight in convergence JSON output.
if ~isstruct(value)
    return;
end
if isfield(value, 'component_moles')
    value.component_moles_total = sum(value.component_moles, 1);
    value = rmfield(value, 'component_moles');
end
if isfield(value, 'mineral_moles')
    value.mineral_moles_total = sum(value.mineral_moles, 1);
    value = rmfield(value, 'mineral_moles');
end
end

function validateFactories(options)
requiredFactories = {'configFactory', 'stateFactory', 'geometryFactory'};
for iFactory = 1:numel(requiredFactories)
    fieldName = requiredFactories{iFactory};
    if ~isfield(options, fieldName) || isempty(options.(fieldName))
        error('RTSPHEM:Benchmark:MissingDriverCaseFactory', ...
            'options.%s is required.', fieldName);
    end
end
end

function tf = shouldUseLightweightStepHistory(options)
tf = logical(getFieldOrDefault(options, 'lightweightStepHistory', false));
end

function summary = runRtSubcyclesLightweight(driver, totalTimeSeconds, ...
        initialState, initialGeometry, cfg, options)
if ~(isscalar(totalTimeSeconds) && isfinite(totalTimeSeconds) && totalTimeSeconds >= 0)
    error('RTSPHEM:Driver:InvalidTimeStep', ...
        'totalTimeSeconds must be a nonnegative finite scalar.');
end
elapsedSeconds = 0;
acceptedSteps = 0;
rejectedSteps = 0;
maxResidual = 0;
maxDisplacement = 0;
dtSeconds = min(initialRtDt(cfg), max(totalTimeSeconds, eps));
maxIterations = getNestedField(cfg, {'time', 'rt', 'maxSubcycles'}, 100000);
observationTimes = getObservationTimes(options);
nextObservation = 1;
observationRecords = struct([]);
lastResult = [];

while elapsedSeconds < totalTimeSeconds - eps
    if acceptedSteps + rejectedSteps >= maxIterations
        error('RTSPHEM:Driver:MaxSubcyclesExceeded', ...
            'RT subcycle loop exceeded maxSubcycles.');
    end
    remainingSeconds = totalTimeSeconds - elapsedSeconds;
    dtTry = min(dtSeconds, remainingSeconds);
    result = driver.runOneStep(dtTry);
    lastResult = result;

    if result.diagnostics.accepted
        acceptedSteps = acceptedSteps + 1;
        elapsedSeconds = elapsedSeconds + dtTry;
        if isfield(result.diagnostics, 'component_absolute_residual_moles')
            maxResidual = max(maxResidual, ...
                max(abs(result.diagnostics.component_absolute_residual_moles(:)), [], 'omitnan'));
        end
        if isfield(result.geometry_info, 'max_displacement_over_h')
            maxDisplacement = max(maxDisplacement, ...
                result.geometry_info.max_displacement_over_h);
        end
        while nextObservation <= numel(observationTimes) && ...
                elapsedSeconds >= observationTimes(nextObservation) - eps
            record = observationRecordFromState(elapsedSeconds, ...
                acceptedSteps + rejectedSteps, result.state, result.geometry, ...
                initialState, initialGeometry, cfg);
            record.requested_time_s = observationTimes(nextObservation);
            observationRecords = appendObservationRecord(observationRecords, record);
            nextObservation = nextObservation + 1;
        end
        dtSeconds = min(maxRtDt(cfg), max(dtTry, eps));
    else
        rejectedSteps = rejectedSteps + 1;
        if isfield(result.solver_state, 'abort') && result.solver_state.abort
            error('RTSPHEM:Driver:SubcycleRetryAborted', ...
                'RT subcycle retry aborted: %s.', string(result.solver_state.abort_reason));
        end
        dtSeconds = result.solver_state.dt_s;
    end
end

summary = struct();
summary.time_s = elapsedSeconds;
summary.accepted_steps = acceptedSteps;
summary.rejected_steps = rejectedSteps;
summary.step_results = struct([]);
if isempty(lastResult)
    summary.state = initialState;
    summary.geometry = initialGeometry;
    summary.solver_state = struct();
else
    summary.state = lastResult.state;
    summary.geometry = lastResult.geometry;
    summary.solver_state = lastResult.solver_state;
end
summary.lightweight_step_history = true;
summary.max_component_mass_residual_moles = maxResidual;
summary.max_displacement_over_h = maxDisplacement;
summary.observation_records = observationRecords;
summary.observation_times_s = observationTimes;
summary.diagnostic_paths = struct();
end

function value = initialRtDt(cfg)
value = getNestedField(cfg, {'time', 'rt', 'initialDt_s'}, 1);
value = min(value, maxRtDt(cfg));
value = max(value, eps);
end

function value = maxRtDt(cfg)
value = getNestedField(cfg, {'time', 'rt', 'maxDt_s'}, Inf);
value = max(value, eps);
end

function connectivity = buildConnectivity(options, refinementScale, runInfo)
if isfield(options, 'connectivityFactory') && ~isempty(options.connectivityFactory)
    connectivity = options.connectivityFactory(refinementScale, runInfo);
else
    connectivity = struct();
end
end

function summary = enrichSummary(driverSummary, initialState, initialGeometry, ...
        cfg, refinementScale, runInfo, options)
summary = driverSummary;
summary.run_name = string(getFieldOrDefault(runInfo, 'name', "driver_case"));
summary.refinement_scale = refinementScale;
summary.accepted = didCompleteRequestedTime(driverSummary, options);
summary.failure_message = "";
summary.initial_porosity = porosityFromGeometry(initialGeometry);
summary.initial_surface_area_cm2 = geometryScalar(initialGeometry, ...
    'initial_surface_area_cm2', 'interface_area_cm2');
summary.initial_solid_volume_cm3 = geometryScalar(initialGeometry, ...
    'initial_solid_volume_cm3', 'solid_volume_cm3');
summary.initial_mineral_moles = sum(initialState.mineral_moles(:), 'omitnan');
summary.final_mineral_moles = sum(driverSummary.state.mineral_moles(:), 'omitnan');
summary.mineral_dissolved_moles = ...
    summary.initial_mineral_moles - summary.final_mineral_moles;
summary.solid_volume_change_cm3 = solidVolumeChangeFromMineral( ...
    summary.mineral_dissolved_moles, cfg);
summary.final_solid_volume_cm3 = max(summary.initial_solid_volume_cm3 + ...
    summary.solid_volume_change_cm3, 0);
summary.final_surface_area_cm2 = geometryScalar(driverSummary.geometry, ...
    'final_surface_area_cm2', 'interface_area_cm2');
summary.final_porosity = finalPorosityFromMineralChange( ...
    initialGeometry, summary.initial_porosity, summary.mineral_dissolved_moles, cfg);
summary.mean_effective_rate_mol_cm2_s = meanEffectiveRate( ...
    summary.mineral_dissolved_moles, summary.initial_surface_area_cm2, ...
    summary.time_s);
summary.initial_component_moles_total = sum(initialState.component_moles, 1);
summary.final_component_moles_total = sum(driverSummary.state.component_moles, 1);
summary.max_component_mass_residual_moles = maxComponentResidual(driverSummary);
summary.max_displacement_over_h = maxDisplacementOverH(driverSummary);
summary.runtime_manifest = rtm.diagnostics.CreateRuntimeManifest(cfg);
summary = addBenchmarkMetadata(summary, cfg, initialGeometry, options);
summary = addObservationRecords(summary, driverSummary, initialState, ...
    initialGeometry, cfg, options);
end

function summary = rejectedSummary(err, initialState, initialGeometry, cfg, ...
        refinementScale, runInfo, options)
summary = struct();
summary.time_s = initialState.time_s;
summary.accepted_steps = 0;
summary.rejected_steps = 1;
summary.step_results = struct([]);
summary.state = initialState;
summary.geometry = initialGeometry;
summary.solver_state = struct('abort', true, ...
    'abort_reason', string(err.message));
summary.run_name = string(getFieldOrDefault(runInfo, 'name', "driver_case"));
summary.refinement_scale = refinementScale;
summary.accepted = false;
summary.failure_message = string(err.message);
summary.initial_porosity = porosityFromGeometry(initialGeometry);
summary.initial_surface_area_cm2 = geometryScalar(initialGeometry, ...
    'initial_surface_area_cm2', 'interface_area_cm2');
summary.initial_solid_volume_cm3 = geometryScalar(initialGeometry, ...
    'initial_solid_volume_cm3', 'solid_volume_cm3');
summary.initial_mineral_moles = sum(initialState.mineral_moles(:), 'omitnan');
summary.final_mineral_moles = summary.initial_mineral_moles;
summary.mineral_dissolved_moles = 0;
summary.solid_volume_change_cm3 = solidVolumeChangeFromMineral( ...
    summary.mineral_dissolved_moles, cfg);
summary.final_solid_volume_cm3 = max(summary.initial_solid_volume_cm3 + ...
    summary.solid_volume_change_cm3, 0);
summary.final_surface_area_cm2 = geometryScalar(initialGeometry, ...
    'final_surface_area_cm2', 'interface_area_cm2');
summary.final_porosity = finalPorosityFromMineralChange( ...
    initialGeometry, summary.initial_porosity, 0, cfg);
summary.mean_effective_rate_mol_cm2_s = meanEffectiveRate( ...
    summary.mineral_dissolved_moles, summary.initial_surface_area_cm2, ...
    summary.time_s);
summary.initial_component_moles_total = sum(initialState.component_moles, 1);
summary.final_component_moles_total = summary.initial_component_moles_total;
summary.max_component_mass_residual_moles = NaN;
summary.max_displacement_over_h = NaN;
summary.runtime_manifest = rtm.diagnostics.CreateRuntimeManifest(cfg);
summary = addBenchmarkMetadata(summary, cfg, initialGeometry, options);
summary = addObservationRecords(summary, summary, initialState, ...
    initialGeometry, cfg, options);
end

function tf = didCompleteRequestedTime(driverSummary, options)
requestedTime = getScalarOption(options, 'totalTime_s', 1);
abort = logical(getNestedField(driverSummary, {'solver_state', 'abort'}, false));
timeTolerance = max(1e-12, 1e-12 .* max(abs(requestedTime), 1));
tf = ~abort && driverSummary.time_s >= requestedTime - timeTolerance;
end

function value = porosityFromGeometry(geometry)
if isstruct(geometry) && isfield(geometry, 'water_volume_cm3') && ...
        isfield(geometry, 'cell_volume_cm3')
    waterVolume = sum(max(geometry.water_volume_cm3(:), 0), 'omitnan');
    cellVolume = sum(max(geometry.cell_volume_cm3(:), 0), 'omitnan');
elseif isstruct(geometry) && isfield(geometry, 'fluid_fraction')
    value = mean(geometry.fluid_fraction(:), 'omitnan');
    return;
else
    value = NaN;
    return;
end
if cellVolume > 0
    value = waterVolume ./ cellVolume;
else
    value = NaN;
end
end

function value = finalPorosityFromMineralChange(geometry, initialPorosity, dissolvedMoles, cfg)
if isstruct(geometry) && isfield(geometry, 'cell_volume_cm3')
    cellVolume = sum(max(geometry.cell_volume_cm3(:), 0), 'omitnan');
else
    cellVolume = NaN;
end
if ~(isfinite(cellVolume) && cellVolume > 0)
    value = initialPorosity;
    return;
end
molarVolume = getNestedField(cfg, {'geometry', 'molarVolume_cm3_mol'}, 1);
waterVolumeChange = max(dissolvedMoles, 0) .* molarVolume;
value = min(max(initialPorosity + waterVolumeChange ./ cellVolume, 0), 1);
end

function value = solidVolumeChangeFromMineral(dissolvedMoles, cfg)
molarVolume = getNestedField(cfg, {'geometry', 'molarVolume_cm3_mol'}, 1);
if ~(isfinite(dissolvedMoles) && isfinite(molarVolume))
    value = NaN;
else
    value = -max(dissolvedMoles, 0) .* molarVolume;
end
end

function value = maxComponentResidual(driverSummary)
if isfield(driverSummary, 'lightweight_step_history') && ...
        logical(driverSummary.lightweight_step_history) && ...
        isfield(driverSummary, 'max_component_mass_residual_moles')
    value = driverSummary.max_component_mass_residual_moles;
    return;
end
value = 0;
for iStep = 1:numel(driverSummary.step_results)
    diagnostics = driverSummary.step_results(iStep).diagnostics;
    if ~(isfield(diagnostics, 'accepted') && logical(diagnostics.accepted))
        continue;
    end
    if isfield(diagnostics, 'component_absolute_residual_moles')
        value = max(value, max(abs(diagnostics.component_absolute_residual_moles(:)), [], 'omitnan'));
    end
end
end

function value = maxDisplacementOverH(driverSummary)
if isfield(driverSummary, 'lightweight_step_history') && ...
        logical(driverSummary.lightweight_step_history) && ...
        isfield(driverSummary, 'max_displacement_over_h')
    value = driverSummary.max_displacement_over_h;
    return;
end
value = 0;
for iStep = 1:numel(driverSummary.step_results)
    diagnostics = driverSummary.step_results(iStep).diagnostics;
    if ~(isfield(diagnostics, 'accepted') && logical(diagnostics.accepted))
        continue;
    end
    geometryInfo = driverSummary.step_results(iStep).geometry_info;
    if isfield(geometryInfo, 'max_displacement_over_h')
        value = max(value, geometryInfo.max_displacement_over_h);
    end
end
end

function summary = addObservationRecords(summary, driverSummary, initialState, ...
        initialGeometry, cfg, options)
observationTimes = getObservationTimes(options);
summary.observation_times_s = observationTimes;
summary.observation_records = struct([]);
if isfield(driverSummary, 'lightweight_step_history') && ...
        logical(driverSummary.lightweight_step_history)
    if isfield(driverSummary, 'observation_records') && ...
            ~isempty(driverSummary.observation_records)
        summary.observation_records = driverSummary.observation_records;
    end
    return;
end
if isempty(observationTimes)
    return;
end

acceptedRecords = acceptedStepRecords(driverSummary, initialState, ...
    initialGeometry, cfg);
if isempty(acceptedRecords)
    return;
end
records = struct([]);
for iObs = 1:numel(observationTimes)
    targetTime = observationTimes(iObs);
    idx = find([acceptedRecords.time_s] >= targetTime - eps, 1);
    if isempty(idx)
        continue;
    end
    record = acceptedRecords(idx);
    record.requested_time_s = targetTime;
    records = appendObservationRecord(records, record);
end
summary.observation_records = records;
end

function observationTimes = getObservationTimes(options)
if isfield(options, 'observationTimes_s') && ~isempty(options.observationTimes_s)
    observationTimes = options.observationTimes_s(:);
else
    observationTimes = zeros(0, 1);
end
if any(~isfinite(observationTimes)) || any(observationTimes < 0)
    error('RTSPHEM:Benchmark:InvalidDriverCaseOption', ...
        'options.observationTimes_s must contain nonnegative finite values.');
end
observationTimes = sort(observationTimes);
end

function records = acceptedStepRecords(driverSummary, initialState, ...
        initialGeometry, cfg)
records = struct([]);
elapsedSeconds = 0;
for iStep = 1:numel(driverSummary.step_results)
    result = driverSummary.step_results(iStep);
    if ~(isfield(result, 'diagnostics') && result.diagnostics.accepted)
        continue;
    end
    dtSeconds = getNestedField(result, {'solver_state', 'dt_s'}, NaN);
    if ~(isscalar(dtSeconds) && isfinite(dtSeconds) && dtSeconds >= 0)
        dtSeconds = 0;
    end
    elapsedSeconds = elapsedSeconds + dtSeconds;
    record = observationRecordFromState(elapsedSeconds, iStep, ...
        result.state, result.geometry, initialState, initialGeometry, cfg);
    records = appendObservationRecord(records, record);
end
end

function record = observationRecordFromState(timeSeconds, stepIndex, state, ...
        geometry, initialState, initialGeometry, cfg)
record = emptyObservationRecord();
record.requested_time_s = NaN;
record.time_s = timeSeconds;
record.source_step_index = stepIndex;
initialMineralMoles = sum(initialState.mineral_moles(:), 'omitnan');
finalMineralMoles = sum(state.mineral_moles(:), 'omitnan');
record.final_mineral_moles = finalMineralMoles;
record.mineral_dissolved_moles = initialMineralMoles - finalMineralMoles;
record.initial_solid_volume_cm3 = geometryScalar(initialGeometry, ...
    'initial_solid_volume_cm3', 'solid_volume_cm3');
record.solid_volume_change_cm3 = solidVolumeChangeFromMineral( ...
    record.mineral_dissolved_moles, cfg);
record.final_solid_volume_cm3 = max(record.initial_solid_volume_cm3 + ...
    record.solid_volume_change_cm3, 0);
record.initial_surface_area_cm2 = geometryScalar(initialGeometry, ...
    'initial_surface_area_cm2', 'interface_area_cm2');
record.final_surface_area_cm2 = geometryScalar(geometry, ...
    'final_surface_area_cm2', 'interface_area_cm2');
record.mean_effective_rate_mol_cm2_s = meanEffectiveRate( ...
    record.mineral_dissolved_moles, record.initial_surface_area_cm2, ...
    timeSeconds);
end

function record = emptyObservationRecord()
record = struct();
record.requested_time_s = NaN;
record.time_s = NaN;
record.source_step_index = NaN;
record.final_mineral_moles = NaN;
record.mineral_dissolved_moles = NaN;
record.initial_solid_volume_cm3 = NaN;
record.solid_volume_change_cm3 = NaN;
record.final_solid_volume_cm3 = NaN;
record.initial_surface_area_cm2 = NaN;
record.final_surface_area_cm2 = NaN;
record.mean_effective_rate_mol_cm2_s = NaN;
end

function records = appendObservationRecord(records, record)
if isempty(records)
    records = record;
else
    records(end + 1, 1) = record;
end
end

function value = meanEffectiveRate(dissolvedMoles, surfaceAreaCm2, elapsedSeconds)
if isfinite(dissolvedMoles) && isfinite(surfaceAreaCm2) && ...
        isfinite(elapsedSeconds) && surfaceAreaCm2 > 0 && elapsedSeconds > 0
    value = dissolvedMoles ./ surfaceAreaCm2 ./ elapsedSeconds;
else
    value = NaN;
end
end

function value = geometryScalar(geometry, preferredField, fallbackField)
if isstruct(geometry) && isfield(geometry, preferredField) && ...
        ~isempty(geometry.(preferredField))
    value = geometry.(preferredField);
    if isscalar(value) && isfinite(value)
        return;
    end
end
if isstruct(geometry) && isfield(geometry, fallbackField) && ...
        ~isempty(geometry.(fallbackField))
    value = sum(geometry.(fallbackField)(:), 'omitnan');
else
    value = NaN;
end
end

function summary = addBenchmarkMetadata(summary, cfg, geometry, options)
mesh = getNestedField(cfg, {'benchmark', 'mesh'}, []);
if isstruct(mesh) && ~isempty(mesh)
    mesh.cell_count = mesh.nx .* mesh.ny;
    summary.benchmark_mesh = mesh;
elseif isstruct(geometry) && isfield(geometry, 'mesh_resolution')
    resolution = geometry.mesh_resolution;
    mesh = struct();
    mesh.nx = resolution(1);
    mesh.ny = resolution(2);
    mesh.cell_count = mesh.nx .* mesh.ny;
    mesh.resolution_label = string(getFieldOrDefault(geometry, ...
        'mesh_resolution_label', sprintf('%dx%d', mesh.nx, mesh.ny)));
    summary.benchmark_mesh = mesh;
end
if isfield(options, 'acceptanceMatrix') && ~isempty(options.acceptanceMatrix)
    summary.acceptance_matrix = options.acceptanceMatrix;
end
end

function value = getScalarOption(options, fieldName, defaultValue)
if isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
if ~(isscalar(value) && isfinite(value) && value >= 0)
    error('RTSPHEM:Benchmark:InvalidDriverCaseOption', ...
        'options.%s must be a nonnegative finite scalar.', fieldName);
end
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end

function value = getNestedField(structValue, path, defaultValue)
value = structValue;
for iPath = 1:numel(path)
    if ~isstruct(value) || ~isfield(value, path{iPath}) || isempty(value.(path{iPath}))
        value = defaultValue;
        return;
    end
    value = value.(path{iPath});
end
end

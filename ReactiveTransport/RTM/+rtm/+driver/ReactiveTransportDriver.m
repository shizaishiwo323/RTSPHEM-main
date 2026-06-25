classdef ReactiveTransportDriver < handle
    %REACTIVETRANSPORTDRIVER Minimal conservative_v2 RTM driver.
    %
    % The initial implementation wires strict_molins reaction-only steps
    % through transaction, chemistry result application, geometry diagnostics,
    % accepted-step validation, and dt retry. Full transport/flow/geometric
    % remap loops are added as separate verified capabilities.

    properties (Access = private)
        Config
        State
        Geometry
        Connectivity
        SolverState
        PhreeqcSession = []
    end

    methods
        function obj = ReactiveTransportDriver(config, state, geometry, connectivity)
            if nargin < 4 || isempty(connectivity)
                connectivity = struct();
            end
            obj.Config = rtm.config.ValidateReactiveTransportConfig(config);
            rtm.state.ValidateState(state);
            obj.State = state;
            obj.Geometry = geometry;
            obj.Connectivity = connectivity;
            obj.SolverState = struct('dt_s', NaN, 'rejected_steps', 0, ...
                'rejection_log', []);
            obj.PhreeqcSession = obj.createPhreeqcSessionIfNeeded();
        end

        function delete(obj)
            if ~isempty(obj.PhreeqcSession)
                try
                    obj.PhreeqcSession.close();
                catch
                end
            end
        end

        function result = runOneStep(obj, dtSeconds)
            if ~(isscalar(dtSeconds) && isfinite(dtSeconds) && dtSeconds > 0)
                error('RTSPHEM:Driver:InvalidTimeStep', ...
                    'dtSeconds must be a positive finite scalar.');
            end
            obj.SolverState.dt_s = dtSeconds;
            tx = rtm.driver.StepTransaction(obj.State, obj.Geometry, obj.SolverState);
            preStepState = obj.State;
            transportLedger = emptyTransportLedger(preStepState, dtSeconds);
            reactionResult = emptyReactionResult(preStepState);
            chemistryLedger = emptyChemistryLedger(preStepState);
            geometryInfo = emptyGeometryInfo(obj.Geometry);

            try
                [transportedState, transportLedger] = obj.transport(dtSeconds);
                obj.State = transportedState;
                reactionInput = obj.buildReactionInput();
                reactionResult = obj.react(reactionInput, dtSeconds);
                [candidateState, chemistryLedger] = rtm.chemistry.ApplyReactionResult( ...
                    obj.State, reactionResult);
                geometryInfo = rtm.geometry.AdvanceGeometryFromMineralMoles( ...
                    obj.Geometry, reactionResult.realized_interface_moles, ...
                    obj.geometryOptions());
                geometryInfo = obj.annotateActualGeometryChange(geometryInfo);
                stepInfo = obj.buildStepInfo(preStepState, candidateState, ...
                    transportLedger, chemistryLedger, reactionResult, geometryInfo);
                diagnostics = rtm.diagnostics.ValidateAcceptedStep( ...
                    stepInfo, obj.diagnosticOptions());

                if diagnostics.accepted
                    obj.State = candidateState;
                    obj.Geometry = rtm.geometry.ApplyMineralVolumeChange( ...
                        obj.Geometry, geometryInfo);
                    obj.SolverState.dt_s = dtSeconds;
                    [stateOut, geometryOut, solverStateOut] = tx.commit( ...
                        obj.State, obj.Geometry, obj.SolverState);
                else
                    [stateOut, geometryOut, solverStateOut] = tx.rollback();
                    solverStateOut = rtm.driver.RejectAndShrinkTimeStep( ...
                        solverStateOut, join(diagnostics.reasons, "; "), obj.retryOptions());
                    obj.State = stateOut;
                    obj.Geometry = geometryOut;
                    obj.SolverState = solverStateOut;
                end
            catch ME
                diagnostics = failedStepDiagnostics(ME);
                [stateOut, geometryOut, solverStateOut] = tx.rollback();
                solverStateOut = rtm.driver.RejectAndShrinkTimeStep( ...
                    solverStateOut, join(diagnostics.reasons, "; "), obj.retryOptions());
                obj.State = stateOut;
                obj.Geometry = geometryOut;
                obj.SolverState = solverStateOut;
            end

            result = struct();
            result.state = stateOut;
            result.geometry = geometryOut;
            result.solver_state = solverStateOut;
            result.reaction_result = reactionResult;
            result.transport_ledger = transportLedger;
            result.chemistry_ledger = chemistryLedger;
            result.geometry_info = geometryInfo;
            result.diagnostics = diagnostics;
            result.transaction_status = tx.status();
        end

        function summary = runRtSubcycles(obj, totalTimeSeconds)
            if ~(isscalar(totalTimeSeconds) && isfinite(totalTimeSeconds) && totalTimeSeconds >= 0)
                error('RTSPHEM:Driver:InvalidTimeStep', ...
                    'totalTimeSeconds must be a nonnegative finite scalar.');
            end
            elapsedSeconds = 0;
            acceptedSteps = 0;
            rejectedSteps = 0;
            stepResults = struct([]);
            dtSeconds = min(obj.initialRtDt(), max(totalTimeSeconds, eps));
            maxIterations = getNestedField(obj.Config, {'time', 'rt', 'maxSubcycles'}, 100000);

            while elapsedSeconds < totalTimeSeconds - eps
                if numel(stepResults) >= maxIterations
                    error('RTSPHEM:Driver:MaxSubcyclesExceeded', ...
                        'RT subcycle loop exceeded maxSubcycles.');
                end
                remainingSeconds = totalTimeSeconds - elapsedSeconds;
                dtTry = min(dtSeconds, remainingSeconds);
                result = obj.runOneStep(dtTry);
                stepResults = appendStepResult(stepResults, result);

                if result.diagnostics.accepted
                    acceptedSteps = acceptedSteps + 1;
                    elapsedSeconds = elapsedSeconds + dtTry;
                    dtSeconds = min(obj.maxRtDt(), max(dtTry, eps));
                else
                    rejectedSteps = rejectedSteps + 1;
                    if isfield(result.solver_state, 'abort') && result.solver_state.abort
                        obj.writeDiagnostics(stepResults);
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
            summary.step_results = stepResults;
            summary.state = obj.State;
            summary.geometry = obj.Geometry;
            summary.solver_state = obj.SolverState;
            summary.diagnostic_paths = obj.writeDiagnostics(stepResults);
        end
    end

    methods (Access = private)
        function reactionInput = buildReactionInput(obj)
            switch obj.Config.chemistry.mode
                case 'strict_molins'
                    options = struct();
                    options.rate_constant_cm_s = getNestedField(obj.Config, ...
                        {'chemistry', 'rate_constant_cm_s'}, 0);
                    options.maxReactantFraction = getNestedField(obj.Config, ...
                        {'time', 'rt', 'maxReactantFraction'}, Inf);
                    options.maxMineralFraction = getNestedField(obj.Config, ...
                        {'time', 'geometry', 'maxMineralFraction'}, Inf);
                    options.reaction_time_integration = getNestedField(obj.Config, ...
                        {'chemistry', 'reaction_time_integration'}, 'explicit_euler');
                    options.reactionClusterDepthCells = getNestedField(obj.Config, ...
                        {'chemistry', 'reaction_cluster_depth_cells'}, 1);
                    reactionInput = rtm.chemistry.BuildStrictMolinsReactionInput( ...
                        obj.State, obj.Geometry, obj.Connectivity, options);
                case {'external_tst_phreeqc', 'phreeqc_kinetics'}
                    reactionInput = getNestedField(obj.Config, ...
                        {'chemistry', 'options'}, struct());
                    if ~isstruct(reactionInput)
                        error('RTSPHEM:Driver:InvalidChemistryOptions', ...
                            'config.chemistry.options must be a struct.');
                    end
                    if ~isfield(reactionInput, 'reactionClusters') || ...
                            isempty(reactionInput.reactionClusters)
                        clusterOptions = struct();
                        clusterOptions.reactionClusterDepthCells = getNestedField( ...
                            obj.Config, {'chemistry', 'reaction_cluster_depth_cells'}, 1);
                        clusterOptions.useReactionClusters = getNestedField( ...
                            obj.Config, {'chemistry', 'useReactionClusters'}, true);
                        reactionInput.reactionClusters = ...
                            rtm.chemistry.BuildReactionClusters(obj.State, ...
                            obj.Geometry, obj.Connectivity, clusterOptions);
                    end
                    if strcmp(obj.Config.chemistry.mode, 'external_tst_phreeqc')
                        reactionInput.rate_constant_cm_s = getNestedField(obj.Config, ...
                            {'chemistry', 'rate_constant_cm_s'}, ...
                            getFieldOrDefault(reactionInput, 'rate_constant_cm_s', 0));
                        reactionInput.maxMineralFraction = getNestedField(obj.Config, ...
                            {'time', 'geometry', 'maxMineralFraction'}, ...
                            getFieldOrDefault(reactionInput, 'maxMineralFraction', Inf));
                    end
                    if ~isempty(obj.PhreeqcSession)
                        reactionInput.phreeqcSession = obj.PhreeqcSession;
                        reactionInput.databasePath = getNestedField(obj.Config, ...
                            {'phreeqc', 'databasePath'}, ...
                            getFieldOrDefault(reactionInput, 'databasePath', ''));
                    end
                otherwise
                    error('RTSPHEM:Driver:UnsupportedChemistryMode', ...
                        ['ReactiveTransportDriver currently supports strict_molins, ', ...
                        'external_tst_phreeqc, and phreeqc_kinetics only.']);
            end
        end

        function [transportedState, transportLedger] = transport(obj, dtSeconds)
            options = getNestedField(obj.Config, {'transport', 'options'}, struct());
            switch obj.Config.transport.backend
                case 'cut_cell_fv'
                    options = normalizeTransportOptions(options, obj.State);
                    [transportedState, transportLedger] = rtm.transport.CutCellTransportFV( ...
                        obj.State, obj.Geometry, dtSeconds, options);
                otherwise
                    error('RTSPHEM:Driver:UnsupportedTransportBackend', ...
                        'ReactiveTransportDriver currently supports cut_cell_fv only.');
            end
        end

        function reactionResult = react(obj, reactionInput, dtSeconds)
            switch obj.Config.chemistry.mode
                case 'strict_molins'
                    reactionResult = rtm.chemistry.StrictMolinsBackend( ...
                        obj.State, obj.Geometry, dtSeconds, reactionInput);
                case 'external_tst_phreeqc'
                    reactionResult = rtm.chemistry.ExternalTstPhreeqcBackend( ...
                        obj.State, obj.Geometry, dtSeconds, reactionInput);
                case 'phreeqc_kinetics'
                    reactionResult = rtm.chemistry.PhreeqcKineticsBackend( ...
                        obj.State, obj.Geometry, dtSeconds, reactionInput);
                otherwise
                    error('RTSPHEM:Driver:UnsupportedChemistryMode', ...
                        ['ReactiveTransportDriver currently supports strict_molins, ', ...
                        'external_tst_phreeqc, and phreeqc_kinetics only.']);
            end
        end

        function stepInfo = buildStepInfo(obj, preStepState, candidateState, ...
                transportLedger, chemistryLedger, reactionResult, geometryInfo)
            oldComponentTotals = sum(preStepState.component_moles, 1);
            newComponentTotals = sum(candidateState.component_moles, 1);
            stepInfo = struct();
            stepInfo.mass = struct();
            stepInfo.mass.initial_component_moles_total = oldComponentTotals;
            stepInfo.mass.final_component_moles_total = newComponentTotals;
            stepInfo.mass.component_residual_moles = ...
                newComponentTotals - oldComponentTotals - ...
                transportLedger.source_delta_moles_total - ...
                transportLedger.internal_flux_delta_moles_total - ...
                getFieldOrDefault(transportLedger, ...
                    'boundary_delta_moles_total', zeros(size(oldComponentTotals))) - ...
                chemistryLedger.component_delta_moles_total;

            stepInfo.geometry = geometryInfo;
            stepInfo.chemistry = struct();
            stepInfo.chemistry.converged = getFieldOrDefault(reactionResult, ...
                'converged', true);
            stepInfo.chemistry.failed_cells = getFieldOrDefault(reactionResult, ...
                'failed_cells', []);
            stepInfo.chemistry.error_message = getFieldOrDefault(reactionResult, ...
                'error_message', "");
            stepInfo.chemistry.charge_balance_residual_eq = getFieldOrDefault( ...
                reactionResult, 'charge_balance_residual_eq', []);
        end

        function options = geometryOptions(obj)
            options = struct();
            options.molarVolume_cm3_mol = getNestedField(obj.Config, ...
                {'geometry', 'molarVolume_cm3_mol'}, 1);
            options.maxDisplacementOverH = getNestedField(obj.Config, ...
                {'geometry', 'maxDisplacementOverH'}, 0.25);
        end

        function options = diagnosticOptions(obj)
            options = struct();
            options.mass_absolute_tolerance_mol = getNestedField(obj.Config, ...
                {'mass', 'absoluteTolerance_mol'}, 1e-14);
            options.mass_relative_tolerance = getNestedField(obj.Config, ...
                {'mass', 'relativeTolerance'}, 1e-8);
            options.solid_absolute_tolerance_cm3 = getNestedField(obj.Config, ...
                {'geometry', 'solidAbsoluteTolerance_cm3'}, 1e-14);
            options.solid_relative_tolerance = getNestedField(obj.Config, ...
                {'geometry', 'solidRelativeTolerance'}, 1e-8);
            options.charge_absolute_tolerance_eq = getNestedField(obj.Config, ...
                {'chemistry', 'chargeAbsoluteTolerance_eq'}, Inf);
            options.max_displacement_over_h = getNestedField(obj.Config, ...
                {'geometry', 'maxDisplacementOverH'}, 0.25);
        end

        function options = retryOptions(obj)
            options = struct();
            options.shrinkFactor = getNestedField(obj.Config, ...
                {'failure', 'shrinkFactor'}, 0.5);
            options.minDt_s = getNestedField(obj.Config, ...
                {'failure', 'minDt_s'}, 1e-8);
            options.maxRetries = getNestedField(obj.Config, ...
                {'failure', 'maxRetries'}, 12);
        end

        function geometryInfo = annotateActualGeometryChange(obj, geometryInfo)
            if isfield(geometryInfo, 'accepted') && ~logical(geometryInfo.accepted)
                solidVolume = obj.Geometry.solid_volume_cm3(:);
                geometryInfo.cell_actual_solid_volume_change_cm3 = zeros(size(solidVolume));
                geometryInfo.actual_solid_volume_change_cm3 = 0;
                geometryInfo.actual_solid_volume_after_cm3 = sum(solidVolume, 'omitnan');
                return;
            end
            candidateGeometry = rtm.geometry.ApplyMineralVolumeChange( ...
                obj.Geometry, geometryInfo);
            oldSolidVolume = obj.Geometry.solid_volume_cm3(:);
            newSolidVolume = candidateGeometry.solid_volume_cm3(:);
            cellActualChange = newSolidVolume - oldSolidVolume;
            geometryInfo.cell_actual_solid_volume_change_cm3 = cellActualChange;
            geometryInfo.actual_solid_volume_change_cm3 = sum(cellActualChange);
            geometryInfo.actual_solid_volume_after_cm3 = sum(newSolidVolume);
        end

        function value = initialRtDt(obj)
            value = getNestedField(obj.Config, {'time', 'rt', 'initialDt_s'}, 1);
            value = min(value, obj.maxRtDt());
            value = max(value, eps);
        end

        function value = maxRtDt(obj)
            value = getNestedField(obj.Config, {'time', 'rt', 'maxDt_s'}, Inf);
            value = max(value, eps);
        end

        function value = diagnosticsOutputDir(obj)
            value = string(getNestedField(obj.Config, ...
                {'output', 'diagnosticsDir'}, ""));
        end

        function paths = writeDiagnostics(obj, stepResults)
            paths = struct();
            diagnosticsDir = obj.diagnosticsOutputDir();
            if strlength(diagnosticsDir) == 0
                return;
            end
            diagnosticsDir = char(diagnosticsDir);
            if exist(diagnosticsDir, 'dir') ~= 7
                mkdir(diagnosticsDir);
            end
            paths = rtm.diagnostics.WriteStepDiagnosticTables( ...
                diagnosticsDir, stepResults);
            manifest = rtm.diagnostics.CreateRuntimeManifest(obj.Config);
            paths.runtime_manifest = ...
                rtm.diagnostics.WriteRuntimeManifest(diagnosticsDir, manifest);
        end

        function session = createPhreeqcSessionIfNeeded(obj)
            session = [];
            if ~ismember(obj.Config.chemistry.mode, ...
                    {'external_tst_phreeqc', 'phreeqc_kinetics'})
                return;
            end
            if ~logical(getNestedField(obj.Config, {'phreeqc', 'persistSession'}, false)) || ...
                    ~logical(getNestedField(obj.Config, {'phreeqc', 'useRunString'}, false))
                return;
            end
            if strcmp(getNestedField(obj.Config, {'phreeqc', 'engine'}, 'none'), 'none')
                return;
            end
            chemistryOptions = getNestedField(obj.Config, {'chemistry', 'options'}, struct());
            if isstruct(chemistryOptions) && ...
                    isfield(chemistryOptions, 'runBatchFunction') && ...
                    ~isempty(chemistryOptions.runBatchFunction)
                return;
            end

            runtimeConfig = struct();
            runtimeConfig.engineType = getNestedField(obj.Config, ...
                {'phreeqc', 'engine'}, 'iphreeqc_com');
            runtimeConfig.comProgId = getNestedField(obj.Config, ...
                {'phreeqc', 'comProgId'}, 'IPhreeqcCOM.Object');
            engineFactory = getNestedField(obj.Config, ...
                {'phreeqc', 'engineFactory'}, []);
            if ~isempty(engineFactory)
                runtimeConfig.engineFactory = engineFactory;
            end
            session = rtm.phreeqc.PhreeqcSession(runtimeConfig);
        end
    end
end

function stepResults = appendStepResult(stepResults, result)
if isempty(stepResults)
    stepResults = result;
else
    stepResults(end + 1) = result;
end
end

function options = normalizeTransportOptions(options, state)
if ~isstruct(options)
    error('RTSPHEM:Driver:InvalidTransportOptions', ...
        'config.transport.options must be a struct.');
end
if isfield(options, 'component_source_mol_s') && isscalar(options.component_source_mol_s)
    options.component_source_mol_s = repmat(options.component_source_mol_s, ...
        size(state.component_moles));
end
end

function ledger = emptyTransportLedger(state, dtSeconds)
numComponents = numel(state.component_names);
ledger = struct();
ledger.dt_s = dtSeconds;
ledger.component_names = state.component_names;
ledger.initial_moles_total = sum(state.component_moles, 1);
ledger.final_moles_total = sum(state.component_moles, 1);
ledger.source_delta_moles_total = zeros(1, numComponents);
ledger.internal_flux_delta_moles_total = zeros(1, numComponents);
ledger.boundary_delta_moles_total = zeros(1, numComponents);
ledger.boundary_advective_delta_moles_total = zeros(1, numComponents);
ledger.boundary_diffusive_delta_moles_total = zeros(1, numComponents);
ledger.component_residual_moles = zeros(1, numComponents);
ledger.max_abs_component_residual_moles = 0;
end

function ledger = emptyChemistryLedger(state)
numComponents = numel(state.component_names);
numMinerals = numel(state.mineral_names);
ledger = struct();
ledger.component_names = state.component_names;
ledger.mineral_names = state.mineral_names;
ledger.initial_component_moles_total = sum(state.component_moles, 1);
ledger.final_component_moles_total = sum(state.component_moles, 1);
ledger.component_delta_moles_total = zeros(1, numComponents);
ledger.initial_mineral_moles_total = sum(state.mineral_moles, 1);
ledger.final_mineral_moles_total = sum(state.mineral_moles, 1);
ledger.mineral_delta_moles_total = zeros(1, numMinerals);
ledger.realized_interface_moles_total = 0;
ledger.converged = false;
ledger.failed_cells = [];
ledger.error_message = "";
end

function result = emptyReactionResult(state)
result = struct();
result.component_delta_moles = zeros(size(state.component_moles));
result.mineral_delta_moles = zeros(size(state.mineral_moles));
result.realized_interface_moles = zeros(size(state.mineral_moles, 1), 1);
result.candidate_interface_moles = zeros(size(state.mineral_moles, 1), 1);
result.converged = false;
result.failed_cells = [];
result.error_message = "";
result.aux = struct('chemistry_mode', "");
end

function geometryInfo = emptyGeometryInfo(geometry)
solidVolume = getVectorFieldOrDefault(geometry, 'solid_volume_cm3', 0);
numCells = numel(solidVolume);
geometryInfo = struct();
geometryInfo.realized_mineral_moles = zeros(numCells, 1);
geometryInfo.molar_volume_cm3_mol = 1;
geometryInfo.solid_volume_before_cm3 = sum(solidVolume, 'omitnan');
geometryInfo.solid_volume_after_cm3 = geometryInfo.solid_volume_before_cm3;
geometryInfo.expected_solid_volume_change_cm3 = 0;
geometryInfo.actual_solid_volume_change_cm3 = 0;
geometryInfo.cell_solid_volume_change_cm3 = zeros(numCells, 1);
geometryInfo.cell_actual_solid_volume_change_cm3 = zeros(numCells, 1);
geometryInfo.max_displacement_over_h = 0;
geometryInfo.accepted = true;
geometryInfo.reject_reason = "";
end

function diagnostics = failedStepDiagnostics(ME)
diagnostics = struct();
diagnostics.accepted = false;
diagnostics.reasons = "state update failed: " + string(ME.message);
diagnostics.component_absolute_residual_moles = NaN;
diagnostics.component_relative_residual = NaN;
diagnostics.solid_volume_absolute_residual_cm3 = NaN;
diagnostics.solid_volume_relative_residual = NaN;
diagnostics.charge_absolute_residual_eq = NaN;
diagnostics.failed_cells = [];
diagnostics.error_message = string(ME.message);
diagnostics.max_displacement_over_h = NaN;
end

function values = getVectorFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    values = structValue.(fieldName)(:);
else
    values = defaultValue;
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

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end

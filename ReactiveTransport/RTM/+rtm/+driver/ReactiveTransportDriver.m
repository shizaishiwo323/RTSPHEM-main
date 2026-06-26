classdef ReactiveTransportDriver < handle
    %REACTIVETRANSPORTDRIVER Minimal conservative_v2 RTM driver.
    %
    % Conservative_v2 RTM driver for PHREEQC-backed calcite chemistry,
    % transaction handling, geometry diagnostics, accepted-step validation,
    % and dt retry.

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
            remapLedger = emptyRemapLedger(preStepState);
            geometryInfo = emptyGeometryInfo(obj.Geometry);
            flowDiagnostics = emptyFlowDiagnostics(size(preStepState.component_moles, 1));

            try
                [transportedState, transportLedger, flowDiagnostics] = obj.transport(dtSeconds);
                obj.State = transportedState;
                reactionInput = obj.buildReactionInput();
                reactionResult = obj.react(reactionInput, dtSeconds);
                [candidateState, chemistryLedger] = rtm.chemistry.ApplyReactionResult( ...
                    obj.State, reactionResult);
                [geometryInfo, candidateGeometry] = obj.advanceGeometryFromMoles( ...
                    obj.Geometry, reactionResult.realized_interface_moles, ...
                    dtSeconds);
                [remappedState, remapLedger] = obj.remapStateAfterGeometryMove( ...
                    candidateState, obj.Geometry, candidateGeometry);
                stepInfo = obj.buildStepInfo(preStepState, remappedState, ...
                    transportLedger, chemistryLedger, remapLedger, ...
                    reactionResult, geometryInfo, flowDiagnostics);
                diagnostics = rtm.diagnostics.ValidateAcceptedStep( ...
                    stepInfo, obj.diagnosticOptions());

                if diagnostics.accepted
                    obj.State = remappedState;
                    obj.Geometry = candidateGeometry;
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
            result.remap_ledger = remapLedger;
            result.geometry_info = geometryInfo;
            result.flow_diagnostics = flowDiagnostics;
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

        function summary = runFixedGeometryRtSubcycles(obj, totalTimeSeconds, options)
            if nargin < 3 || isempty(options)
                options = struct();
            end
            if ~(isscalar(totalTimeSeconds) && isfinite(totalTimeSeconds) && totalTimeSeconds >= 0)
                error('RTSPHEM:Driver:InvalidTimeStep', ...
                    'totalTimeSeconds must be a nonnegative finite scalar.');
            end
            freezeGeometry = logical(getFieldOrDefault(options, ...
                'freezeGeometry', true));
            freezeMineralInventory = logical(getFieldOrDefault(options, ...
                'freezeMineralInventory', true));

            fixedGeometry = obj.Geometry;
            fixedMineralMoles = obj.State.mineral_moles;
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
                result = obj.runOneStepForFixedGeometry(dtTry, freezeGeometry);

                if result.diagnostics.accepted
                    result = obj.applyFixedGeometryContract(result, ...
                        fixedGeometry, fixedMineralMoles, freezeGeometry, ...
                        freezeMineralInventory);
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
                stepResults = appendStepResult(stepResults, result);
            end

            summary = struct();
            summary.time_s = elapsedSeconds;
            summary.accepted_steps = acceptedSteps;
            summary.rejected_steps = rejectedSteps;
            summary.step_results = stepResults;
            summary.state = obj.State;
            summary.geometry = obj.Geometry;
            summary.solver_state = obj.SolverState;
            summary.fixed_geometry = freezeGeometry;
            summary.fixed_mineral_inventory = freezeMineralInventory;
            summary.diagnostic_paths = obj.writeDiagnostics(stepResults);
        end

        function summary = runGeometryMacroStep(obj, dtGeometrySeconds, options)
            if nargin < 3 || isempty(options)
                options = struct();
            end
            if ~(isscalar(dtGeometrySeconds) && isfinite(dtGeometrySeconds) && ...
                    dtGeometrySeconds >= 0)
                error('RTSPHEM:Driver:InvalidTimeStep', ...
                    'dtGeometrySeconds must be a nonnegative finite scalar.');
            end

            geometryBefore = obj.Geometry;
            tx = rtm.driver.StepTransaction(obj.State, obj.Geometry, ...
                obj.SolverState);
            try
                rtOptions = options;
                rtOptions.freezeGeometry = true;
                rtOptions.freezeMineralInventory = logical(getFieldOrDefault( ...
                    options, 'freezeMineralInventoryDuringRt', false));
                rtSummary = obj.runFixedGeometryRtSubcycles(dtGeometrySeconds, rtOptions);
                accumulatedMoles = accumulatedRealizedMoles(rtSummary, ...
                    numel(getVectorFieldOrDefault(geometryBefore, ...
                    'solid_volume_cm3', 0)));
                [geometryInfo, candidateGeometry] = obj.advanceGeometryFromMoles( ...
                    geometryBefore, accumulatedMoles, dtGeometrySeconds);
                if ~geometryInfo.accepted
                    [stateOut, geometryOut, solverStateOut] = tx.rollback();
                    solverStateOut.dt_s = dtGeometrySeconds;
                    solverStateOut = rtm.driver.RejectAndShrinkTimeStep( ...
                        solverStateOut, geometryInfo.reject_reason, obj.retryOptions());
                    obj.State = stateOut;
                    obj.Geometry = geometryOut;
                    obj.SolverState = solverStateOut;
                    error('RTSPHEM:Driver:GeometryMacroStepRejected', ...
                        'Geometry macro step rejected: %s.', string(geometryInfo.reject_reason));
                end
                [remappedState, remapLedger] = obj.remapStateAfterGeometryMove( ...
                    obj.State, geometryBefore, candidateGeometry);
                obj.State = remappedState;
                obj.Geometry = candidateGeometry;
                obj.SolverState.dt_s = dtGeometrySeconds;
                [~, ~, ~] = tx.commit(obj.State, obj.Geometry, obj.SolverState);
            catch ME
                if ~strcmp(tx.status(), 'rolled_back') && ~strcmp(tx.status(), 'committed')
                    [stateOut, geometryOut, solverStateOut] = tx.rollback();
                    obj.State = stateOut;
                    obj.Geometry = geometryOut;
                    obj.SolverState = solverStateOut;
                end
                rethrow(ME);
            end

            summary = rtSummary;
            summary.rt_summary = rtSummary;
            summary.geometry_before = geometryBefore;
            summary.geometry = obj.Geometry;
            summary.state = obj.State;
            summary.remap_ledger = remapLedger;
            summary.geometry_info = geometryInfo;
            summary.accumulated_realized_mineral_moles = accumulatedMoles;
            summary.geometry_macro_step_s = dtGeometrySeconds;
        end

        function summary = runQuasiSteadyGeometry(obj, totalTimeSeconds, options)
            if nargin < 3 || isempty(options)
                options = struct();
            end
            if ~(isscalar(totalTimeSeconds) && isfinite(totalTimeSeconds) && ...
                    totalTimeSeconds >= 0)
                error('RTSPHEM:Driver:InvalidTimeStep', ...
                    'totalTimeSeconds must be a nonnegative finite scalar.');
            end

            elapsedSeconds = 0;
            acceptedMacroSteps = 0;
            rejectedMacroSteps = 0;
            macroStepResults = struct([]);
            macroStepSizes = zeros(0, 1);
            maxIterations = getNestedField(obj.Config, ...
                {'time', 'geometry', 'maxMacroSteps'}, 100000);
            dtSeconds = min(obj.maxGeometryDt(), max(totalTimeSeconds, eps));

            while elapsedSeconds < totalTimeSeconds - eps
                if acceptedMacroSteps + rejectedMacroSteps >= maxIterations
                    error('RTSPHEM:Driver:MaxGeometryMacroStepsExceeded', ...
                        'Geometry macro loop exceeded maxMacroSteps.');
                end
                remainingSeconds = totalTimeSeconds - elapsedSeconds;
                dtTry = min(dtSeconds, remainingSeconds);
                try
                    stepSummary = obj.runGeometryMacroStep(dtTry, options);
                    macroStepResults = appendStepResult(macroStepResults, ...
                        stepSummary);
                    macroStepSizes(end + 1, 1) = dtTry; %#ok<AGROW>
                    acceptedMacroSteps = acceptedMacroSteps + 1;
                    elapsedSeconds = elapsedSeconds + dtTry;
                    dtSeconds = min(obj.maxGeometryDt(), ...
                        max(totalTimeSeconds - elapsedSeconds, eps));
                catch ME
                    if ~strcmp(ME.identifier, ...
                            'RTSPHEM:Driver:GeometryMacroStepRejected')
                        rethrow(ME);
                    end
                    rejectedMacroSteps = rejectedMacroSteps + 1;
                    if isfield(obj.SolverState, 'abort') && obj.SolverState.abort
                        error('RTSPHEM:Driver:GeometryMacroRetryAborted', ...
                            'Geometry macro retry aborted: %s.', ...
                            string(obj.SolverState.abort_reason));
                    end
                    dtSeconds = obj.currentRetryDt();
                end
            end

            summary = struct();
            summary.time_s = elapsedSeconds;
            summary.accepted_macro_steps = acceptedMacroSteps;
            summary.rejected_macro_steps = rejectedMacroSteps;
            summary.macro_step_results = macroStepResults;
            summary.macro_step_sizes_s = macroStepSizes;
            summary.state = obj.State;
            summary.geometry = obj.Geometry;
            summary.solver_state = obj.SolverState;
        end
    end

    methods (Access = private)
        function reactionInput = buildReactionInput(obj)
            switch obj.Config.chemistry.mode
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
                    reactionInput.calciteStoichiometryAbsoluteTolerance_mol = ...
                        getNestedField(obj.Config, ...
                        {'chemistry', 'calciteStoichiometryAbsoluteTolerance_mol'}, ...
                        getFieldOrDefault(reactionInput, ...
                        'calciteStoichiometryAbsoluteTolerance_mol', 1e-14));
                    reactionInput.calciteStoichiometryRelativeTolerance = ...
                        getNestedField(obj.Config, ...
                        {'chemistry', 'calciteStoichiometryRelativeTolerance'}, ...
                        getFieldOrDefault(reactionInput, ...
                        'calciteStoichiometryRelativeTolerance', 1e-8));
                    if ~isempty(obj.PhreeqcSession)
                        reactionInput.phreeqcSession = obj.PhreeqcSession;
                        reactionInput.databasePath = getNestedField(obj.Config, ...
                            {'phreeqc', 'databasePath'}, ...
                            getFieldOrDefault(reactionInput, 'databasePath', ''));
                    end
                otherwise
                    error('RTSPHEM:Driver:UnsupportedChemistryMode', ...
                        ['ReactiveTransportDriver currently supports ', ...
                        'external_tst_phreeqc and phreeqc_kinetics only.']);
            end
        end

        function [transportedState, transportLedger, flowDiagnostics] = transport(obj, dtSeconds)
            options = getNestedField(obj.Config, {'transport', 'options'}, struct());
            flow = obj.resolveFlowFaceFluxes();
            flowDiagnostics = obj.flowDiagnostics(flow);
            options = rtm.flow.ApplyFaceFluxesToTransportOptions( ...
                options, flow);
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

        function flow = resolveFlowFaceFluxes(obj)
            flow = getNestedField(obj.Config, {'flow', 'face_fluxes'}, []);
            hyphmOptions = getNestedField(obj.Config, {'flow', 'hyphmStokes'}, []);
            if isstruct(hyphmOptions) && ~isempty(hyphmOptions)
                flow = rtm.flow.SolveHyphmStokesLevel(obj.Geometry, hyphmOptions);
            end
        end

        function diagnostics = flowDiagnostics(obj, flow)
            numCells = size(obj.State.component_moles, 1);
            if ~isstruct(flow) || isempty(fieldnames(flow))
                diagnostics = emptyFlowDiagnostics(numCells);
                return;
            end
            options = struct();
            options.active_fluid_cell = getVectorFieldOrDefault( ...
                obj.Geometry, 'active_fluid_cell', true(numCells, 1));
            options.absoluteTolerance_cm3_s = getNestedField(obj.Config, ...
                {'flow', 'absoluteTolerance_cm3_s'}, 1e-12);
            options.relativeTolerance = getNestedField(obj.Config, ...
                {'flow', 'relativeTolerance'}, Inf);
            diagnostics = rtm.flow.ValidateFaceFluxDivergence(flow, numCells, options);
        end

        function reactionResult = react(obj, reactionInput, dtSeconds)
            switch obj.Config.chemistry.mode
                case 'external_tst_phreeqc'
                    reactionResult = rtm.chemistry.ExternalTstPhreeqcBackend( ...
                        obj.State, obj.Geometry, dtSeconds, reactionInput);
                case 'phreeqc_kinetics'
                    reactionResult = rtm.chemistry.PhreeqcKineticsBackend( ...
                        obj.State, obj.Geometry, dtSeconds, reactionInput);
                otherwise
                    error('RTSPHEM:Driver:UnsupportedChemistryMode', ...
                        ['ReactiveTransportDriver currently supports ', ...
                        'external_tst_phreeqc and phreeqc_kinetics only.']);
            end
        end

        function stepInfo = buildStepInfo(obj, preStepState, candidateState, ...
                transportLedger, chemistryLedger, remapLedger, reactionResult, ...
                geometryInfo, flowDiagnostics)
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
            stepInfo.transport = transportLedger;
            stepInfo.flow = flowDiagnostics;
            stepInfo.reaction = chemistryLedger;
            stepInfo.remap = struct('component_delta_moles_total', ...
                getFieldOrDefault(remapLedger, 'component_residual_moles', ...
                zeros(size(oldComponentTotals))));

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

        function [geometryInfo, candidateGeometry] = advanceGeometryFromMoles( ...
                obj, geometryBefore, realizedMoles, dtSeconds)
            mesh = obj.geometryMesh();
            if ~isempty(mesh) && isfield(geometryBefore, 'level_set') && ...
                    ~isempty(geometryBefore.level_set)
                [~, geometryInfo] = rtm.geometry.AdvanceLevelSetFromMineralMoles( ...
                    mesh, geometryBefore.level_set, geometryBefore, ...
                    realizedMoles, max(dtSeconds, eps), obj.geometryOptions());
                if geometryInfo.accepted
                    candidateGeometry = geometryInfo.new_geometry;
                else
                    candidateGeometry = geometryBefore;
                    geometryInfo = annotateRejectedActualChange(geometryInfo, ...
                        geometryBefore);
                end
                return;
            end

            geometryInfo = rtm.geometry.AdvanceGeometryFromMineralMoles( ...
                geometryBefore, realizedMoles, obj.geometryOptions());
            if geometryInfo.accepted
                geometryInfo = annotateActualGeometryChangeForGeometry( ...
                    geometryBefore, geometryInfo);
                candidateGeometry = rtm.geometry.ApplyMineralVolumeChange( ...
                    geometryBefore, geometryInfo);
            else
                geometryInfo = annotateRejectedActualChange(geometryInfo, ...
                    geometryBefore);
                candidateGeometry = geometryBefore;
            end
        end

        function [remappedState, remapLedger] = remapStateAfterGeometryMove( ...
                obj, stateBeforeRemap, oldGeometry, newGeometry)
            remapOptions = struct();
            remapOptions.new_mineral_moles = stateBeforeRemap.mineral_moles;
            configuredOptions = getNestedField(obj.Config, ...
                {'geometry', 'remapOptions'}, struct());
            if isstruct(configuredOptions) && isfield(configuredOptions, ...
                    'overlap_volume_cm3')
                remapOptions.overlap_volume_cm3 = configuredOptions.overlap_volume_cm3;
            end
            [remappedState, remapLedger] = rtm.geometry.ConservativeRemap( ...
                stateBeforeRemap, oldGeometry, newGeometry, remapOptions);
        end

        function mesh = geometryMesh(obj)
            mesh = [];
            if isstruct(obj.Connectivity) && isfield(obj.Connectivity, 'mesh') && ...
                    ~isempty(obj.Connectivity.mesh)
                mesh = obj.Connectivity.mesh;
                return;
            end
            mesh = getNestedField(obj.Config, {'geometry', 'mesh'}, []);
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

        function value = maxGeometryDt(obj)
            value = getNestedField(obj.Config, {'time', 'geometry', 'maxDt_s'}, Inf);
            value = max(value, eps);
        end

        function value = currentRetryDt(obj)
            value = getFieldOrDefault(obj.SolverState, 'dt_s', obj.maxGeometryDt());
            if ~(isscalar(value) && isfinite(value) && value > 0)
                value = obj.maxGeometryDt();
            end
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

        function result = applyFixedGeometryContract(obj, result, fixedGeometry, ...
                fixedMineralMoles, freezeGeometry, freezeMineralInventory)
            if freezeMineralInventory
                obj.State.mineral_moles = fixedMineralMoles;
                result.state.mineral_moles = fixedMineralMoles;
                if isfield(result, 'chemistry_ledger')
                    result.chemistry_ledger.final_mineral_moles_total = ...
                        sum(fixedMineralMoles, 1);
                    result.chemistry_ledger.mineral_delta_moles_total = ...
                        zeros(size(result.chemistry_ledger.mineral_delta_moles_total));
                end
            end
            if freezeGeometry
                obj.Geometry = fixedGeometry;
                result.geometry = fixedGeometry;
                result.geometry_info = freezeGeometryInfo(result.geometry_info, ...
                    fixedGeometry);
            end
        end

        function result = runOneStepForFixedGeometry(obj, dtSeconds, freezeGeometry)
            if ~freezeGeometry
                result = obj.runOneStep(dtSeconds);
                return;
            end
            originalConfig = obj.Config;
            restoreConfig = onCleanup(@() obj.restoreConfig(originalConfig));
            if ~isfield(obj.Config, 'geometry') || ~isstruct(obj.Config.geometry)
                obj.Config.geometry = struct();
            end
            obj.Config.geometry.maxDisplacementOverH = Inf;
            result = obj.runOneStep(dtSeconds);
            clear restoreConfig;
            obj.restoreConfig(originalConfig);
        end

        function restoreConfig(obj, config)
            obj.Config = config;
        end
    end
end

function geometryInfo = freezeGeometryInfo(geometryInfo, fixedGeometry)
solidVolume = getVectorFieldOrDefault(fixedGeometry, 'solid_volume_cm3', 0);
geometryInfo.cell_actual_solid_volume_change_cm3 = zeros(size(solidVolume));
geometryInfo.actual_solid_volume_change_cm3 = 0;
geometryInfo.actual_solid_volume_after_cm3 = sum(solidVolume, 'omitnan');
geometryInfo.solid_volume_after_cm3 = sum(solidVolume, 'omitnan');
geometryInfo.fixed_geometry = true;
end

function moles = accumulatedRealizedMoles(summary, numCells)
moles = zeros(numCells, 1);
if ~isfield(summary, 'step_results')
    return;
end
for iStep = 1:numel(summary.step_results)
    result = summary.step_results(iStep);
    if ~(isfield(result, 'diagnostics') && result.diagnostics.accepted)
        continue;
    end
    if ~isfield(result, 'reaction_result') || ...
            ~isfield(result.reaction_result, 'realized_interface_moles')
        continue;
    end
    stepMoles = result.reaction_result.realized_interface_moles(:);
    if numel(stepMoles) ~= numCells
        error('RTSPHEM:Driver:GeometryMacroStepSizeMismatch', ...
            'realized_interface_moles must match the geometry cell count.');
    end
    moles = moles + stepMoles;
end
end

function geometryInfo = annotateActualGeometryChangeForGeometry(geometry, geometryInfo)
candidateGeometry = rtm.geometry.ApplyMineralVolumeChange(geometry, geometryInfo);
oldSolidVolume = geometry.solid_volume_cm3(:);
newSolidVolume = candidateGeometry.solid_volume_cm3(:);
cellActualChange = newSolidVolume - oldSolidVolume;
geometryInfo.cell_actual_solid_volume_change_cm3 = cellActualChange;
geometryInfo.actual_solid_volume_change_cm3 = sum(cellActualChange);
geometryInfo.actual_solid_volume_after_cm3 = sum(newSolidVolume);
end

function geometryInfo = annotateRejectedActualChange(geometryInfo, geometry)
solidVolume = getVectorFieldOrDefault(geometry, 'solid_volume_cm3', 0);
geometryInfo.cell_actual_solid_volume_change_cm3 = zeros(size(solidVolume));
geometryInfo.actual_solid_volume_change_cm3 = 0;
geometryInfo.actual_solid_volume_after_cm3 = sum(solidVolume, 'omitnan');
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

function ledger = emptyRemapLedger(state)
numComponents = numel(state.component_names);
ledger = struct();
ledger.initial_component_moles_total = sum(state.component_moles, 1);
ledger.final_component_moles_total = sum(state.component_moles, 1);
ledger.assigned_component_moles = sum(state.component_moles, 1);
ledger.unassigned_component_moles = zeros(1, numComponents);
ledger.component_residual_moles = zeros(1, numComponents);
ledger.max_abs_component_residual_moles = 0;
ledger.new_water_without_overlap_cells = [];
ledger.old_water_unmapped_cells = [];
ledger.overlap_volume_cm3 = sparse(size(state.component_moles, 1), ...
    size(state.component_moles, 1));
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

function diagnostics = emptyFlowDiagnostics(numCells)
diagnostics = struct();
diagnostics.accepted = true;
diagnostics.cell_divergence_cm3_s = zeros(numCells, 1);
diagnostics.max_abs_cell_divergence_cm3_s = 0;
diagnostics.global_residual_cm3_s = 0;
diagnostics.inlet_flow_cm3_s = 0;
diagnostics.outlet_flow_cm3_s = 0;
diagnostics.inlet_outlet_relative_residual = 0;
diagnostics.absolute_tolerance_cm3_s = Inf;
diagnostics.relative_tolerance = Inf;
diagnostics.active_fluid_cell = true(numCells, 1);
diagnostics.failure_reasons = strings(0, 1);
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

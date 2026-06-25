function paths = WriteStepDiagnosticTables(outputDir, stepResults)
%WRITESTEPDIAGNOSTICTABLES Write conservative RTM step diagnostics to CSV.

outputDir = char(outputDir);
if exist(outputDir, 'dir') ~= 7
    error('RTSPHEM:Diagnostics:MissingOutputDirectory', ...
        'Output directory does not exist: %s.', outputDir);
end
if nargin < 2 || isempty(stepResults)
    stepResults = struct([]);
end

massPath = fullfile(outputDir, 'mass_balance_components.csv');
solidPath = fullfile(outputDir, 'solid_geometry_balance.csv');
chemistryPath = fullfile(outputDir, 'chemistry_step_status.csv');
rejectionPath = fullfile(outputDir, 'step_rejection_log.csv');
fluxPath = fullfile(outputDir, 'transport_flux_balance.csv');
clusterPath = fullfile(outputDir, 'reaction_cluster_diagnostics.csv');

massTable = buildMassTable(stepResults);
solidTable = buildSolidTable(stepResults);
chemistryTable = buildChemistryTable(stepResults);
rejectionTable = buildRejectionTable(stepResults);
fluxTable = buildTransportFluxTable(stepResults);
clusterTable = buildReactionClusterTable(stepResults);

if isempty(massTable), massTable = emptyMassTable(); end
if isempty(solidTable), solidTable = emptySolidTable(); end
if isempty(chemistryTable), chemistryTable = emptyChemistryTable(); end
if isempty(rejectionTable), rejectionTable = emptyRejectionTable(); end
if isempty(fluxTable), fluxTable = emptyTransportFluxTable(); end
if isempty(clusterTable), clusterTable = emptyReactionClusterTable(); end

writetable(massTable, massPath);
writetable(solidTable, solidPath);
writetable(chemistryTable, chemistryPath);
writetable(rejectionTable, rejectionPath);
writetable(fluxTable, fluxPath);
writetable(clusterTable, clusterPath);

paths = struct();
paths.mass_balance_components = string(massPath);
paths.solid_geometry_balance = string(solidPath);
paths.chemistry_step_status = string(chemistryPath);
paths.step_rejection_log = string(rejectionPath);
paths.transport_flux_balance = string(fluxPath);
paths.reaction_cluster_diagnostics = string(clusterPath);
end

function tableOut = buildMassTable(stepResults)
rows = {};
for iStep = 1:numel(stepResults)
    result = stepResults(iStep);
    ledger = requireStructField(result, 'transport_ledger');
    fullLedger = getNestedField(result, {'diagnostics', 'full_mass_ledger'}, struct());
    if isFullMassLedger(fullLedger)
        componentNames = string(getFieldOrDefault(fullLedger, ...
            'component_names', ledger.component_names));
    else
        componentNames = string(ledger.component_names(:));
    end
    componentNames = componentNames(:);
    numComponents = numel(componentNames);
    if isFullMassLedger(fullLedger)
        initialMoles = vectorField(fullLedger, ...
            'initial_component_moles_total', numComponents, 0);
        finalMoles = vectorField(fullLedger, ...
            'final_component_moles_total', numComponents, 0);
        sourceDelta = vectorField(fullLedger, ...
            'transport_source_delta_moles_total', numComponents, 0);
        internalDelta = vectorField(fullLedger, ...
            'transport_internal_flux_delta_moles_total', numComponents, 0);
        boundaryDelta = vectorField(fullLedger, ...
            'transport_boundary_delta_moles_total', numComponents, 0);
        chemistryDelta = vectorField(fullLedger, ...
            'reaction_delta_moles_total', numComponents, 0);
        remapDelta = vectorField(fullLedger, ...
            'remap_delta_moles_total', numComponents, 0);
        residual = vectorField(fullLedger, ...
            'component_residual_moles', numComponents, 0);
    else
        initialMoles = vectorField(ledger, 'initial_moles_total', numComponents, 0);
        finalMoles = vectorField(ledger, 'final_moles_total', numComponents, 0);
        sourceDelta = vectorField(ledger, 'source_delta_moles_total', numComponents, 0);
        internalDelta = vectorField(ledger, 'internal_flux_delta_moles_total', numComponents, 0);
        boundaryDelta = vectorField(ledger, 'boundary_delta_moles_total', numComponents, 0);
        chemistryDelta = vectorField(getStructField(result, 'chemistry_ledger'), ...
            'component_delta_moles_total', numComponents, 0);
        remapDelta = zeros(numComponents, 1);
        residual = massResidual(result, initialMoles, finalMoles, sourceDelta, ...
            internalDelta, boundaryDelta, chemistryDelta, numComponents);
    end
    accepted = logical(getNestedField(result, {'diagnostics', 'accepted'}, false));
    dtSeconds = scalarField(ledger, 'dt_s', NaN);

    for iComponent = 1:numComponents
        rows(end + 1, :) = {iStep, accepted, dtSeconds, ...
            componentNames(iComponent), initialMoles(iComponent), ...
            finalMoles(iComponent), sourceDelta(iComponent), ...
            internalDelta(iComponent), boundaryDelta(iComponent), ...
            chemistryDelta(iComponent), remapDelta(iComponent), ...
            residual(iComponent)}; %#ok<AGROW>
    end
end
if isempty(rows)
    tableOut = table();
else
    tableOut = cell2table(rows, 'VariableNames', massVariableNames());
end
end

function tableOut = buildSolidTable(stepResults)
rows = {};
for iStep = 1:numel(stepResults)
    result = stepResults(iStep);
    geometry = requireStructField(result, 'geometry_info');
    fullLedger = getNestedField(result, {'diagnostics', 'full_mass_ledger'}, struct());
    accepted = logical(getNestedField(result, {'diagnostics', 'accepted'}, false));
    solidBefore = scalarField(geometry, 'solid_volume_before_cm3', NaN);
    solidAfter = scalarField(geometry, 'solid_volume_after_cm3', NaN);
    expectedChange = scalarField(geometry, 'expected_solid_volume_change_cm3', NaN);
    actualChange = scalarField(geometry, 'actual_solid_volume_change_cm3', ...
        expectedChange);
    residual = actualChange - expectedChange;
    realizedMoles = sum(vectorField(geometry, 'realized_mineral_moles', [], 0), ...
        'omitnan');
    initialMinerals = sumVectorFieldOrNaN(fullLedger, ...
        'initial_mineral_moles_total');
    reactionMineralDelta = sumVectorFieldOrNaN(fullLedger, ...
        'reaction_mineral_delta_moles_total');
    finalMinerals = sumVectorFieldOrNaN(fullLedger, ...
        'final_mineral_moles_total');
    mineralResidual = sumVectorFieldOrNaN(fullLedger, ...
        'mineral_residual_moles');
    displacementOverH = scalarField(geometry, 'max_displacement_over_h', NaN);
    rows(end + 1, :) = {iStep, accepted, solidBefore, solidAfter, ...
        expectedChange, actualChange, residual, realizedMoles, ...
        initialMinerals, reactionMineralDelta, finalMinerals, ...
        mineralResidual, displacementOverH}; %#ok<AGROW>
end
if isempty(rows)
    tableOut = table();
else
    tableOut = cell2table(rows, 'VariableNames', solidVariableNames());
end
end

function tableOut = buildChemistryTable(stepResults)
rows = {};
for iStep = 1:numel(stepResults)
    result = stepResults(iStep);
    reaction = getStructField(result, 'reaction_result');
    diagnostics = getStructField(result, 'diagnostics');
    aux = getStructField(reaction, 'aux');
    failedCells = vectorField(reaction, 'failed_cells', [], []);
    rows(end + 1, :) = {iStep, ...
        logical(getFieldOrDefault(diagnostics, 'accepted', false)), ...
        logical(getFieldOrDefault(reaction, 'converged', true)), ...
        numel(failedCells), join(string(failedCells(:)), ";"), ...
        string(getFieldOrDefault(reaction, 'error_message', "")), ...
        scalarField(diagnostics, 'charge_absolute_residual_eq', ...
            chemistryChargeResidual(reaction)), ...
        string(getFieldOrDefault(aux, 'chemistry_mode', "")), ...
        string(getFieldOrDefault(aux, 'phreeqc_run_method', "")), ...
        logical(getFieldOrDefault(aux, 'phreeqc_session_reused', false)), ...
        scalarField(aux, 'phreeqc_run_status', NaN), ...
        string(getFieldOrDefault(result, 'transaction_status', ""))}; %#ok<AGROW>
end
if isempty(rows)
    tableOut = table();
else
    tableOut = cell2table(rows, 'VariableNames', chemistryVariableNames());
end
end

function tableOut = buildRejectionTable(stepResults)
rows = {};
for iStep = 1:numel(stepResults)
    result = stepResults(iStep);
    diagnostics = getStructField(result, 'diagnostics');
    if logical(getFieldOrDefault(diagnostics, 'accepted', false))
        continue;
    end
    ledger = getStructField(result, 'transport_ledger');
    reasons = string(getFieldOrDefault(diagnostics, 'reasons', ""));
    reasons = reasons(strlength(reasons) > 0);
    if isempty(reasons)
        reasons = string(getFieldOrDefault(result, 'transaction_status', "rejected"));
    end
    rows(end + 1, :) = {iStep, scalarField(ledger, 'dt_s', NaN), ...
        join(reasons(:), "; "), ...
        string(getFieldOrDefault(result, 'transaction_status', "")), ...
        scalarField(diagnostics, 'component_absolute_residual_moles', NaN), ...
        scalarField(diagnostics, 'solid_volume_absolute_residual_cm3', NaN), ...
        scalarField(diagnostics, 'mineral_absolute_residual_moles', NaN), ...
        scalarField(diagnostics, 'charge_absolute_residual_eq', NaN)}; %#ok<AGROW>
end
if isempty(rows)
    tableOut = table();
else
    tableOut = cell2table(rows, 'VariableNames', rejectionVariableNames());
end
end

function tableOut = buildTransportFluxTable(stepResults)
rows = {};
for iStep = 1:numel(stepResults)
    result = stepResults(iStep);
    ledger = requireStructField(result, 'transport_ledger');
    componentNames = string(ledger.component_names(:));
    numComponents = numel(componentNames);
    sourceDelta = vectorField(ledger, 'source_delta_moles_total', numComponents, 0);
    internalDelta = vectorField(ledger, 'internal_flux_delta_moles_total', numComponents, 0);
    boundaryDelta = vectorField(ledger, 'boundary_delta_moles_total', numComponents, 0);
    boundaryAdvectiveDelta = vectorField(ledger, ...
        'boundary_advective_delta_moles_total', numComponents, 0);
    boundaryDiffusiveDelta = vectorField(ledger, ...
        'boundary_diffusive_delta_moles_total', numComponents, 0);
    roundoffSuppressed = roundoffSuppressedByComponent(ledger, numComponents);
    residual = vectorField(ledger, 'component_residual_moles', numComponents, 0);
    flowDiagnostics = getStructField(result, 'flow_diagnostics');
    inletFlow = scalarField(flowDiagnostics, 'inlet_flow_cm3_s', NaN);
    outletFlow = scalarField(flowDiagnostics, 'outlet_flow_cm3_s', NaN);
    inletOutletRelativeResidual = scalarField(flowDiagnostics, ...
        'inlet_outlet_relative_residual', NaN);
    accepted = logical(getNestedField(result, {'diagnostics', 'accepted'}, false));
    dtSeconds = scalarField(ledger, 'dt_s', NaN);
    for iComponent = 1:numComponents
        rows(end + 1, :) = {iStep, accepted, dtSeconds, ...
            componentNames(iComponent), sourceDelta(iComponent), ...
            internalDelta(iComponent), boundaryDelta(iComponent), ...
            boundaryAdvectiveDelta(iComponent), boundaryDiffusiveDelta(iComponent), ...
            roundoffSuppressed(iComponent), residual(iComponent), ...
            inletFlow, outletFlow, inletOutletRelativeResidual}; %#ok<AGROW>
    end
end
if isempty(rows)
    tableOut = table();
else
    tableOut = cell2table(rows, 'VariableNames', transportFluxVariableNames());
end
end

function tableOut = buildReactionClusterTable(stepResults)
rows = {};
for iStep = 1:numel(stepResults)
    result = stepResults(iStep);
    reaction = getStructField(result, 'reaction_result');
    diagnostics = getStructField(result, 'diagnostics');
    realized = vectorField(reaction, 'realized_interface_moles', [], 0);
    candidate = vectorField(reaction, 'candidate_interface_moles', [], 0);
    failedCells = vectorField(reaction, 'failed_cells', [], []);
    clusterIds = vectorField(reaction, 'cluster_ids', [], []);
    clusterCount = numel(unique(clusterIds(~isnan(clusterIds))));
    if clusterCount == 0 && any(realized ~= 0 | candidate ~= 0)
        clusterCount = 1;
    end
    rows(end + 1, :) = {iStep, ...
        logical(getFieldOrDefault(diagnostics, 'accepted', false)), ...
        clusterCount, sum(candidate, 'omitnan'), sum(realized, 'omitnan'), ...
        numel(failedCells), string(getFieldOrDefault(reaction, ...
        'error_message', ""))}; %#ok<AGROW>
end
if isempty(rows)
    tableOut = table();
else
    tableOut = cell2table(rows, 'VariableNames', reactionClusterVariableNames());
end
end

function values = massResidual(result, initialMoles, finalMoles, sourceDelta, ...
        internalDelta, boundaryDelta, chemistryDelta, numComponents)
ledger = getStructField(result, 'transport_ledger');
if isfield(ledger, 'component_residual_moles')
    values = vectorField(ledger, 'component_residual_moles', numComponents, 0);
    return;
end
diagnosticResidual = getNestedField(result, ...
    {'diagnostics', 'component_absolute_residual_moles'}, []);
if ~isempty(diagnosticResidual)
    values = vectorWithLength(diagnosticResidual, numComponents, 0);
    return;
end
values = finalMoles - initialMoles - sourceDelta - internalDelta - ...
    boundaryDelta - chemistryDelta;
end

function tableOut = emptyMassTable()
tableOut = table('Size', [0, numel(massVariableNames())], ...
    'VariableTypes', {'double', 'logical', 'double', 'string', ...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
    'double'}, ...
    'VariableNames', massVariableNames());
end

function tableOut = emptySolidTable()
tableOut = table('Size', [0, numel(solidVariableNames())], ...
    'VariableTypes', {'double', 'logical', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double'}, ...
    'VariableNames', solidVariableNames());
end

function tableOut = emptyChemistryTable()
tableOut = table('Size', [0, numel(chemistryVariableNames())], ...
    'VariableTypes', {'double', 'logical', 'logical', 'double', ...
    'string', 'string', 'double', 'string', 'string', 'logical', ...
    'double', 'string'}, ...
    'VariableNames', chemistryVariableNames());
end

function tableOut = emptyRejectionTable()
tableOut = table('Size', [0, numel(rejectionVariableNames())], ...
    'VariableTypes', {'double', 'double', 'string', 'string', 'double', ...
    'double', 'double', 'double'}, ...
    'VariableNames', rejectionVariableNames());
end

function tableOut = emptyTransportFluxTable()
tableOut = table('Size', [0, numel(transportFluxVariableNames())], ...
    'VariableTypes', {'double', 'logical', 'double', 'string', ...
    'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double'}, ...
    'VariableNames', transportFluxVariableNames());
end

function tableOut = emptyReactionClusterTable()
tableOut = table('Size', [0, numel(reactionClusterVariableNames())], ...
    'VariableTypes', {'double', 'logical', 'double', 'double', ...
    'double', 'double', 'string'}, ...
    'VariableNames', reactionClusterVariableNames());
end

function names = massVariableNames()
names = {'step_index', 'accepted', 'dt_s', 'component_name', ...
    'initial_moles', 'final_moles', 'source_delta_moles', ...
    'internal_flux_delta_moles', 'boundary_delta_moles', ...
    'chemistry_delta_moles', 'remap_delta_moles', ...
    'component_residual_moles'};
end

function names = solidVariableNames()
names = {'step_index', 'accepted', 'solid_volume_before_cm3', ...
    'solid_volume_after_cm3', 'expected_solid_volume_change_cm3', ...
    'actual_solid_volume_change_cm3', 'solid_volume_residual_cm3', ...
    'realized_mineral_moles', 'initial_mineral_moles', ...
    'reaction_mineral_delta_moles', 'final_mineral_moles', ...
    'mineral_residual_moles', 'max_displacement_over_h'};
end

function names = chemistryVariableNames()
names = {'step_index', 'accepted', 'converged', 'failed_cell_count', ...
    'failed_cells', 'error_message', 'charge_balance_max_abs_eq', 'chemistry_mode', ...
    'phreeqc_run_method', 'phreeqc_session_reused', 'phreeqc_run_status', ...
    'transaction_status'};
end

function names = rejectionVariableNames()
names = {'step_index', 'dt_s', 'reason', 'transaction_status', ...
    'component_absolute_residual_moles', ...
    'solid_volume_absolute_residual_cm3', ...
    'mineral_absolute_residual_moles', ...
    'charge_absolute_residual_eq'};
end

function names = transportFluxVariableNames()
names = {'step_index', 'accepted', 'dt_s', 'component_name', ...
    'source_delta_moles', 'internal_flux_delta_moles', ...
    'boundary_delta_moles', 'boundary_advective_delta_moles', ...
    'boundary_diffusive_delta_moles', ...
    'roundoff_suppressed_moles', ...
    'component_residual_moles', ...
    'inlet_flow_cm3_s', 'outlet_flow_cm3_s', ...
    'inlet_outlet_relative_residual'};
end

function names = reactionClusterVariableNames()
names = {'step_index', 'accepted', 'cluster_count', ...
    'candidate_interface_moles_total', 'realized_interface_moles_total', ...
    'failed_cell_count', 'error_message'};
end

function values = roundoffSuppressedByComponent(ledger, numComponents)
if isfield(ledger, 'roundoff_suppressed_moles') && ...
        ~isempty(ledger.roundoff_suppressed_moles)
    matrix = ledger.roundoff_suppressed_moles;
    if size(matrix, 2) == numComponents
        values = sum(matrix, 1, 'omitnan').';
        return;
    end
end
values = vectorField(ledger, 'roundoff_suppressed_moles_total', ...
    numComponents, 0);
end

function value = requireStructField(structValue, fieldName)
if ~isstruct(structValue) || ~isfield(structValue, fieldName) || ...
        ~isstruct(structValue.(fieldName))
    error('RTSPHEM:Diagnostics:InvalidStepResult', ...
        'stepResult.%s is required.', fieldName);
end
value = structValue.(fieldName);
end

function value = getStructField(structValue, fieldName)
if isstruct(structValue) && isfield(structValue, fieldName) && ...
        isstruct(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = struct();
end
end

function value = isFullMassLedger(value)
value = isstruct(value) && isfield(value, 'initial_component_moles_total') && ...
    isfield(value, 'final_component_moles_total');
end

function value = getFieldOrDefault(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ~isempty(structValue.(fieldName))
    value = structValue.(fieldName);
else
    value = defaultValue;
end
end

function values = vectorField(structValue, fieldName, expectedLength, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ...
        ~isempty(structValue.(fieldName))
    values = structValue.(fieldName);
else
    values = defaultValue;
end
values = vectorWithLength(values, expectedLength, defaultValue);
end

function values = vectorWithLength(values, expectedLength, defaultValue)
values = values(:);
if isempty(expectedLength)
    return;
end
if isscalar(values) && expectedLength > 1
    values = repmat(values, expectedLength, 1);
elseif numel(values) ~= expectedLength
    values = repmat(defaultValue, expectedLength, 1);
end
end

function value = scalarField(structValue, fieldName, defaultValue)
if isstruct(structValue) && isfield(structValue, fieldName) && ...
        ~isempty(structValue.(fieldName))
    values = structValue.(fieldName);
    value = values(1);
else
    value = defaultValue;
end
end

function value = sumVectorFieldOrNaN(structValue, fieldName)
if ~isstruct(structValue) || ~isfield(structValue, fieldName) || ...
        isempty(structValue.(fieldName))
    value = NaN;
    return;
end
values = structValue.(fieldName);
values = values(:);
values = values(isfinite(values));
if isempty(values)
    value = NaN;
else
    value = sum(values);
end
end

function value = chemistryChargeResidual(reaction)
if ~isstruct(reaction) || ~isfield(reaction, 'charge_balance_residual_eq') || ...
        isempty(reaction.charge_balance_residual_eq)
    value = 0;
    return;
end
values = reaction.charge_balance_residual_eq(:);
values = values(isfinite(values));
if isempty(values)
    value = 0;
else
    value = max(abs(values));
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

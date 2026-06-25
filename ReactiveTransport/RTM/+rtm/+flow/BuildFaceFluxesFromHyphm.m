function flow = BuildFaceFluxesFromHyphm(hyphmResult, numCells, options)
%BUILDFACEFLUXESFROMHYPHM Build RTM transport face fields from HyPHM output.

if nargin < 3 || isempty(options)
    options = struct();
end
flow = rtm.flow.MapRt0FluxToCutCellFaces(hyphmResult);
if nargin >= 2 && ~isempty(numCells)
    flow.num_cells = numCells;
end

validateDivergence = logical(getOption(options, 'validateDivergence', false));
if validateDivergence
    diagnosticsOptions = struct();
    diagnosticsOptions.absoluteTolerance_cm3_s = getOption(options, ...
        'absoluteTolerance_cm3_s', 1e-12);
    if isfield(options, 'active_fluid_cell') && ~isempty(options.active_fluid_cell)
        diagnosticsOptions.active_fluid_cell = options.active_fluid_cell;
    end
    flow.diagnostics = rtm.flow.ValidateFaceFluxDivergence( ...
        flow, numCells, diagnosticsOptions);
    if isfield(options, 'errorOnDivergence') && logical(options.errorOnDivergence) && ...
            ~flow.diagnostics.accepted
        error('RTSPHEM:Flow:FaceFluxDivergence', ...
            'HyPHM face flux divergence exceeds tolerance.');
    end
else
    flow.diagnostics = struct();
end
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

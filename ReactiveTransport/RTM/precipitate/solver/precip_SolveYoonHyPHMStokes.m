function flow = precip_SolveYoonHyPHMStokes(state, spec, options)
% precip_SolveYoonHyPHMStokes - Bridge Yoon masks to an injected HyPHM Stokes solve.
%
% This function is the production-backend seam for Milestone 4. It is
% fail-closed: without an injected HyPHM solver function it errors instead of
% returning proxy data that could pass quantitative readiness gates.

if nargin < 3
    options = struct();
end

solverFcn = getHyphmSolverFcn(spec, options);
levelSetPackage = precip_MapCellMaskToLevelSet(state, spec, options);
hyphmOutput = solverFcn(levelSetPackage, spec, options);
mapped = precip_MapHyPHMFluxToCartesianFaces(hyphmOutput, state, spec, options);

baselineOutletFlux = computeBaselineOutletFlux(state, spec);
relativePermeability = mapped.outletFlux_cm3_s / ...
    max(baselineOutletFlux, eps);
relativePermeability = min(max(relativePermeability, 0), 1);

flow = struct();
flow.solver = 'hyphm_stokes';
flow.isProxy = false;
flow.isStokes = true;
flow.velocityX_cm_s = mapped.velocityX_cm_s;
flow.velocityY_cm_s = mapped.velocityY_cm_s;
flow.inletFlux_cm3_s = mapped.inletFlux_cm3_s;
flow.outletFlux_cm3_s = mapped.outletFlux_cm3_s;
flow.baselineOutletFlux_cm3_s = baselineOutletFlux;
flow.boundaryFluxClosureRelativeError = ...
    mapped.boundaryFluxClosureRelativeError;
flow.maxDivergenceResidual_s_inv = mapped.maxDivergenceResidual_s_inv;
flow.relativePermeability = relativePermeability;
flow.pressureDropRelative = getFieldOrDefault(hyphmOutput, ...
    'pressureDropRelative', 1 / max(relativePermeability, eps));
flow.flowRateRelative = relativePermeability;
flow.numBlockedCells = mapped.numBlockedCells;
flow.numFluidCells = mapped.numFluidCells;
flow.numMobileCells = mapped.numMobileCells;
flow.linearResidualRelative = getFieldOrDefault(hyphmOutput, ...
    'linearResidualRelative', NaN);
flow.source = getFieldOrDefault(mapped, 'source', 'hyphm_stokes');
flow.mappingMethod = mapped.mappingMethod;
flow.levelSetSignConvention = levelSetPackage.levelSetSignConvention;
flow.note = ['Injected HyPHM Stokes bridge; production acceptance still ', ...
    'requires residual, divergence, and boundary-closure evidence.'];
end

function solverFcn = getHyphmSolverFcn(spec, options)
solverFcn = [];
if isstruct(options) && isfield(options, 'hyphmStokesSolverFcn') && ...
        ~isempty(options.hyphmStokesSolverFcn)
    solverFcn = options.hyphmStokesSolverFcn;
elseif isfield(spec, 'hyphmStokesSolverFcn') && ...
        ~isempty(spec.hyphmStokesSolverFcn)
    solverFcn = spec.hyphmStokesSolverFcn;
end
if isempty(solverFcn)
    error('RTSPHEM:Precipitate:MissingHyphmStokesSolver', ...
        ['spec.hyphmStokesSolverFcn or options.hyphmStokesSolverFcn ', ...
        'is required for yoonFlowSolver = hyphm_stokes.']);
end
if ~isa(solverFcn, 'function_handle')
    error('RTSPHEM:Precipitate:InvalidHyphmStokesSolver', ...
        'HyPHM Stokes solver must be a function handle.');
end
end

function flux = computeBaselineOutletFlux(state, spec)
substrateMask = logicalMask(state, 'substrateMask');
dy = getFieldOrDefault(spec, 'dy_cm', 1);
depth = getFieldOrDefault(spec, 'thickness_cm', ...
    getFieldOrDefault(spec, 'depth_cm', 1));
velocity = getFieldOrDefault(spec, 'darcyVelocity_cm_s', 0);
flux = velocity * dy * depth * nnz(~substrateMask(:, end));
end

function mask = logicalMask(state, fieldName)
if isfield(state, fieldName) && ~isempty(state.(fieldName))
    mask = logical(state.(fieldName));
elseif isfield(state, 'Vm')
    mask = false(size(state.Vm));
else
    error('RTSPHEM:Precipitate:MissingVm', ...
        'state.Vm is required when %s is absent.', fieldName);
end
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

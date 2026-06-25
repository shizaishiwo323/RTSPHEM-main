function flow = precip_RecomputeYoonFlowField(state, spec, flowSolverFcn)
% precip_RecomputeYoonFlowField - Recompute or proxy Yoon flow diagnostics.
%
% Inputs:
%   state         - Yoon state with substrateMask and blockedMask.
%   spec          - benchmark spec passed through to external solvers.
%   flowSolverFcn - optional function handle: flow = f(state, spec).
%
% Output:
%   flow - struct with relative permeability and flow redistribution fields.

if nargin < 3
    flowSolverFcn = [];
end

if ~isempty(flowSolverFcn)
    if ~isa(flowSolverFcn, 'function_handle')
        error('RTSPHEM:Precipitate:InvalidFlowSolverFcn', ...
            'flowSolverFcn must be a function handle.');
    end
    flow = flowSolverFcn(state, spec);
    flow = normalizeFlowStruct(flow, state, false);
    return;
end

solverName = lower(strrep(char(string(getFieldOrDefault(spec, ...
    'yoonFlowSolver', 'proxy_blocked_fraction'))), '-', '_'));
switch solverName
    case {'proxy_blocked_fraction', 'proxy'}
    case {'finite_volume_darcy_pressure', 'darcy_pressure', ...
            'fv_darcy_pressure'}
        flow = precip_SolveYoonDarcyFlow2D(state, spec);
        flow = normalizeFlowStruct(flow, state, false);
        return;
    case {'finite_difference_stokes_brinkman', 'stokes_brinkman', ...
            'fd_stokes_brinkman'}
        flow = precip_SolveYoonStokesBrinkman2D(state, spec);
        flow = normalizeFlowStruct(flow, state, false);
        return;
    otherwise
        error('RTSPHEM:Precipitate:InvalidYoonFlowSolver', ...
            'Unsupported Yoon flow solver: %s.', solverName);
end

blockedMask = logical(getStateMask(state, 'blockedMask'));
substrateMask = logical(getStateMask(state, 'substrateMask'));
mobileMask = ~substrateMask;
numFluidCells = nnz(mobileMask);
numBlockedCells = nnz(blockedMask & mobileMask);
blockedFraction = numBlockedCells / max(numFluidCells, 1);
relativePermeability = max(0, 1 - blockedFraction);

flow = struct();
flow.solver = 'proxy_blocked_fraction';
flow.isProxy = true;
flow.isStokes = false;
flow.blockedFraction = blockedFraction;
flow.relativePermeability = relativePermeability;
flow.pressureDropRelative = 1 / max(relativePermeability, eps);
flow.flowRateRelative = relativePermeability;
flow.numBlockedCells = numBlockedCells;
flow.numFluidCells = numFluidCells;
flow.note = ['Proxy diagnostic only; attach flowSolverFcn for Stokes ', ...
    'recomputation.'];
end

function mask = getStateMask(state, fieldName)
if isfield(state, fieldName) && ~isempty(state.(fieldName))
    mask = state.(fieldName);
elseif isfield(state, 'Vm')
    mask = false(size(state.Vm));
else
    error('RTSPHEM:Precipitate:MissingVm', ...
        'state.Vm is required when %s is absent.', fieldName);
end
end

function flow = normalizeFlowStruct(flow, state, isProxy)
if ~isstruct(flow)
    error('RTSPHEM:Precipitate:InvalidFlowSolverOutput', ...
        'flowSolverFcn must return a struct.');
end
blockedMask = logical(getStateMask(state, 'blockedMask'));
substrateMask = logical(getStateMask(state, 'substrateMask'));
mobileMask = ~substrateMask;
numFluidCells = nnz(mobileMask);
numBlockedCells = nnz(blockedMask & mobileMask);
blockedFraction = numBlockedCells / max(numFluidCells, 1);

if ~isfield(flow, 'solver') || isempty(flow.solver)
    flow.solver = 'external_flow_solver';
end
flow.isProxy = isProxy;
flow.isStokes = getFieldOrDefault(flow, 'isStokes', false);
flow.numBlockedCells = getFieldOrDefault(flow, 'numBlockedCells', ...
    numBlockedCells);
flow.numFluidCells = getFieldOrDefault(flow, 'numFluidCells', ...
    numFluidCells);
flow.blockedFraction = getFieldOrDefault(flow, 'blockedFraction', ...
    blockedFraction);
flow.relativePermeability = getFieldOrDefault(flow, ...
    'relativePermeability', max(0, 1 - flow.blockedFraction));
flow.pressureDropRelative = getFieldOrDefault(flow, ...
    'pressureDropRelative', 1 / max(flow.relativePermeability, eps));
flow.flowRateRelative = getFieldOrDefault(flow, 'flowRateRelative', ...
    flow.relativePermeability);
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

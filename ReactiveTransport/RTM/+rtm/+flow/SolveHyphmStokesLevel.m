function flow = SolveHyphmStokesLevel(geometry, options)
%SOLVEHYPHMSTOKESLEVEL Adapter from a HyPHM/Stokes solve to RTM face fluxes.
%
% The actual Stokes discretization remains supplied by solverFunction. This
% adapter standardizes the returned RT0/HyPHM fluxes for conservative_v2.

if nargin < 2 || isempty(options)
    options = struct();
end
if ~isstruct(options) || ~isfield(options, 'solverFunction') || ...
        isempty(options.solverFunction)
    error('RTSPHEM:Flow:MissingHyphmSolver', ...
        'options.solverFunction is required to solve HyPHM Stokes flow.');
end
solverFunction = options.solverFunction;
if ~isa(solverFunction, 'function_handle')
    error('RTSPHEM:Flow:InvalidHyphmSolver', ...
        'options.solverFunction must be a function handle.');
end

hyphmResult = solverFunction(geometry, options);
numCells = inferCellCount(geometry, options);
flow = rtm.flow.BuildFaceFluxesFromHyphm(hyphmResult, numCells, options);
flow.solver_type = "hyphm_stokes";
end

function numCells = inferCellCount(geometry, options)
if isfield(options, 'numCells') && ~isempty(options.numCells)
    numCells = options.numCells;
elseif isstruct(geometry) && isfield(geometry, 'water_volume_cm3') && ...
        ~isempty(geometry.water_volume_cm3)
    numCells = numel(geometry.water_volume_cm3);
else
    error('RTSPHEM:Flow:InvalidGeometry', ...
        'Cannot infer Stokes cell count from geometry.');
end
if ~(isscalar(numCells) && isfinite(numCells) && numCells == round(numCells) && ...
        numCells >= 1)
    error('RTSPHEM:Flow:InvalidGeometry', ...
        'Stokes cell count must be a positive integer scalar.');
end
end

function result = precip_YoonMicrocontinuumSolver(command, spec, options)
% precip_YoonMicrocontinuumSolver - Minimal Yoon Vm micro-continuum entry point.
%
% Inputs:
%   command - currently supports 'initialize'.
%   spec    - struct from precip_ZhangYoonBenchmarkSpec.
%   options - reserved optional struct.
%
% Output:
%   result  - initialized solver state with static substrate and dynamic Vm.

if nargin < 3
    options = struct();
end
if nargin < 2 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end
if ~ischar(command) && ~isstring(command)
    error('RTSPHEM:Precipitate:InvalidYoonSolverCommand', ...
        'Command must be text.');
end

switch lower(char(command))
    case 'initialize'
        result = initializeState(spec, options);
    otherwise
        error('RTSPHEM:Precipitate:InvalidYoonSolverCommand', ...
            'Unsupported Yoon solver command: %s.', char(command));
end
end

function state = initializeState(spec, options)
[xCenters, yCenters] = cellCenters(spec);
[substrateMask, geometry] = buildSubstrateMask(spec, xCenters, yCenters, ...
    options);
vm = zeros(spec.numY, spec.numX);

components = struct();
components.names = spec.componentNames;
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    components.(fieldName) = initialComponentField(spec, fieldName);
end

state = struct();
state.modelFamily = spec.modelFamily;
state.grid.xCenters_cm = xCenters;
state.grid.yCenters_cm = yCenters;
state.grid.dx_cm = spec.dx_cm;
state.grid.dy_cm = spec.dy_cm;
state.substrateMask = substrateMask;
state.substrateGeometrySource = geometry.source;
state.substrateMaskFile = geometry.filePath;
state.Vm = vm;
state.blockedMask = vm >= spec.blockedVmThreshold;
state.components = components;
state.fluidVolumeFraction = double(~substrateMask) .* (1 - vm);
end

function [xCenters, yCenters] = cellCenters(spec)
xCenters = ((1:spec.numX) - 0.5) .* spec.dx_cm;
yCenters = ((1:spec.numY)' - 0.5) .* spec.dy_cm;
end

function [substrateMask, geometry] = buildSubstrateMask(spec, xCenters, ...
    yCenters, options)
if isfield(options, 'substrateMask')
    substrateMask = logical(options.substrateMask);
    if ~isequal(size(substrateMask), [spec.numY, spec.numX])
        error('RTSPHEM:Precipitate:InvalidSubstrateMask', ...
            'substrateMask must match [numY, numX].');
    end
    geometry = struct('source', 'options.substrateMask', 'filePath', '');
    return;
end
if isfield(spec, 'substrateMaskFile') && ~isempty(spec.substrateMaskFile)
    geometry = precip_LoadZhangYoonGeometry(spec.substrateMaskFile, spec);
    substrateMask = geometry.substrateMask;
    return;
end

% Keep the first Yoon path independent from dynamic CaCO3. The static
% substrate mask is present from the start and can later be replaced by a
% digitized Zhang/Yoon geometry without changing Vm semantics.
substrateMask = false(spec.numY, spec.numX);
if ~isfield(spec, 'postDiameter_cm') || isempty(spec.postDiameter_cm)
    geometry = struct('source', 'empty_substrate', 'filePath', '');
    return;
end

[xGrid, yGrid] = meshgrid(xCenters, yCenters);
radius = 0.5 * spec.postDiameter_cm;
pitchX = getFieldOrDefault(spec, 'postPitchX_cm', 0.048);
pitchY = getFieldOrDefault(spec, 'postPitchY_cm', 0.034);
for x0 = pitchX / 2:pitchX:spec.lengthXAxis_cm
    rowOffset = 0.5 * pitchY * mod(round(x0 / pitchX), 2);
    for y0 = pitchY / 2 + rowOffset:pitchY:spec.lengthYAxis_cm
        substrateMask = substrateMask | ((xGrid - x0).^2 + (yGrid - y0).^2 <= radius.^2);
    end
end
geometry = struct('source', 'approximate_cylindrical_post_array', ...
    'filePath', '');
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function values = initialComponentField(spec, fieldName)
values = ones(spec.numY, spec.numX) .* spec.initial.(fieldName);
end

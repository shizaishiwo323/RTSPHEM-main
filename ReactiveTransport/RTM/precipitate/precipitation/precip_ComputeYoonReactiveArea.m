function areaCm2 = precip_ComputeYoonReactiveArea(state, spec)
% precip_ComputeYoonReactiveArea - Yoon reactive area per fluid cell.
%
% Inputs:
%   state - Yoon state with substrateMask and optional Vm.
%   spec  - benchmark spec with dx/dy/thickness.
%
% Output:
%   areaCm2 - per-cell reactive area, cm2.

if ~isfield(state, 'substrateMask')
    error('RTSPHEM:Precipitate:MissingSubstrateMask', ...
        'state.substrateMask is required.');
end

substrateMask = logical(state.substrateMask);
filledMask = filledPrecipitateMask(state, spec);
fluidMask = ~(substrateMask | filledMask);

areaCm2 = zeros(size(substrateMask));
areaCm2(fluidMask) = 2 * spec.dx_cm * spec.dy_cm;
areaCm2 = areaCm2 + exposedVerticalFaceArea(fluidMask, substrateMask, spec);
areaCm2 = areaCm2 + exposedVerticalFaceArea(fluidMask, filledMask, spec);
areaCm2(~fluidMask) = 0;
end

function filledMask = filledPrecipitateMask(state, spec)
if ~isfield(state, 'Vm') || isempty(state.Vm)
    filledMask = false(size(state.substrateMask));
    return;
end
vm = state.Vm;
if ~isequal(size(vm), size(state.substrateMask))
    error('RTSPHEM:Precipitate:InvalidVmSize', ...
        'state.Vm must match state.substrateMask size.');
end
tolerance = getFieldOrDefault(spec, 'filledVmTolerance', 1e-12);
filledMask = vm >= 1 - tolerance;
filledMask(state.substrateMask) = false;
end

function area = exposedVerticalFaceArea(fluidMask, solidMask, spec)
area = zeros(size(fluidMask));
if isempty(fluidMask)
    return;
end

horizontalFaceArea = spec.dy_cm * spec.thickness_cm;
verticalFaceArea = spec.dx_cm * spec.thickness_cm;

leftNeighborSolid = [false(size(solidMask, 1), 1), solidMask(:, 1:end-1)];
rightNeighborSolid = [solidMask(:, 2:end), false(size(solidMask, 1), 1)];
downNeighborSolid = [false(1, size(solidMask, 2)); solidMask(1:end-1, :)];
upNeighborSolid = [solidMask(2:end, :); false(1, size(solidMask, 2))];

area = area + horizontalFaceArea .* fluidMask .* ...
    (double(leftNeighborSolid) + double(rightNeighborSolid));
area = area + verticalFaceArea .* fluidMask .* ...
    (double(downNeighborSolid) + double(upNeighborSolid));
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

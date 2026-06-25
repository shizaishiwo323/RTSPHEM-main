function metrics = precip_ComputeYoonAreaMetrics(state, spec, masks)
% precip_ComputeYoonAreaMetrics - Compute Yoon precipitated-area diagnostics.
%
% Inputs:
%   state - Yoon state with Vm and optional substrateMask.
%   spec  - benchmark spec with dx/dy and areaVmThreshold.
%   masks - optional struct with firstPoreMask and firstThreePoresMask.
%
% Output:
%   metrics - area metrics based on Vm > areaVmThreshold.

if nargin < 3
    masks = struct();
end
if ~isfield(state, 'Vm')
    error('RTSPHEM:Precipitate:MissingVm', 'state.Vm is required.');
end

activeMask = state.Vm > spec.areaVmThreshold;
if isfield(state, 'substrateMask')
    activeMask = activeMask & ~state.substrateMask;
end
cellAreaCm2 = spec.dx_cm * spec.dy_cm;

metrics = struct();
metrics.vmThreshold = spec.areaVmThreshold;
metrics.totalPrecipitatedArea_cm2 = nnz(activeMask) * cellAreaCm2;
metrics.firstPorePrecipitatedArea_cm2 = regionArea(activeMask, masks, ...
    'firstPoreMask', cellAreaCm2);
metrics.firstThreePoresPrecipitatedArea_cm2 = regionArea(activeMask, masks, ...
    'firstThreePoresMask', cellAreaCm2);
metrics.numActiveCells = nnz(activeMask);
end

function area = regionArea(activeMask, masks, fieldName, cellAreaCm2)
if isfield(masks, fieldName) && ~isempty(masks.(fieldName))
    regionMask = logical(masks.(fieldName));
    if ~isequal(size(regionMask), size(activeMask))
        error('RTSPHEM:Precipitate:InvalidYoonAreaMask', ...
            '%s must match Vm size.', fieldName);
    end
    area = nnz(activeMask & regionMask) * cellAreaCm2;
else
    area = NaN;
end
end

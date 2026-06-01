function cfg = GetFlowBoundaryConfig(flowDirection, lengthXAxis, lengthYAxis, nxParts, nyParts)
% GetFlowBoundaryConfig - Boundary and axis mapping for the RTM flow direction.
%
% Boundary IDs in the HyPHM grid are:
%   1 = bottom, 2 = right, 3 = top, 4 = left.

if nargin < 1 || isempty(flowDirection)
    flowDirection = 'left_to_right';
end

direction = lower(strtrim(char(flowDirection)));
direction = strrep(direction, '-', '_');

switch direction
    case {'left_to_right', 'x_positive', '+x'}
        cfg.direction = 'left_to_right';
        cfg.inletId = 4;
        cfg.outletId = 2;
        cfg.wallIds = [1, 3];
        cfg.axisIndex = 1;
        cfg.transverseAxisIndex = 2;
        cfg.flowLength = lengthXAxis;
        cfg.crossSectionLength = lengthYAxis;
        cfg.axisParts = nxParts;
        cfg.velocityVector = @(speed) [speed; 0];
        cfg.inletPredicate = @(x) x(1) < eps;
        cfg.outletEdgePredicate = @(baryE) baryE(:, 1) > (lengthXAxis - eps);
        cfg.inletNodePredicate = @(coordV) coordV(:, 1) < eps;
        cfg.outletNodePredicate = @(coordV) coordV(:, 1) > (lengthXAxis - eps);
        cfg.outletTrianglePredicate = @(baryT) baryT(:, 1) > (lengthXAxis - 2 * lengthXAxis / nxParts);
        cfg.axisCoordinate = @(coord) coord(:, 1);

    case {'bottom_to_top', 'up', 'y_positive', '+y'}
        cfg.direction = 'bottom_to_top';
        cfg.inletId = 1;
        cfg.outletId = 3;
        cfg.wallIds = [4, 2];
        cfg.axisIndex = 2;
        cfg.transverseAxisIndex = 1;
        cfg.flowLength = lengthYAxis;
        cfg.crossSectionLength = lengthXAxis;
        cfg.axisParts = nyParts;
        cfg.velocityVector = @(speed) [0; speed];
        cfg.inletPredicate = @(x) x(2) < eps;
        cfg.outletEdgePredicate = @(baryE) baryE(:, 2) > (lengthYAxis - eps);
        cfg.inletNodePredicate = @(coordV) coordV(:, 2) < eps;
        cfg.outletNodePredicate = @(coordV) coordV(:, 2) > (lengthYAxis - eps);
        cfg.outletTrianglePredicate = @(baryT) baryT(:, 2) > (lengthYAxis - 2 * lengthYAxis / nyParts);
        cfg.axisCoordinate = @(coord) coord(:, 2);

    otherwise
        error('MATLAB:GetFlowBoundaryConfig:UnknownDirection', ...
            'Unknown flowDirection "%s". Use left_to_right or bottom_to_top.', flowDirection);
end

cfg.stokesDirichletIds = [cfg.inletId, cfg.wallIds];
cfg.transportNeumannIds = [cfg.wallIds, cfg.outletId];
cfg.transportFluxIds = cfg.inletId;
cfg.segmentEdges = linspace(0, cfg.flowLength, 6);
end

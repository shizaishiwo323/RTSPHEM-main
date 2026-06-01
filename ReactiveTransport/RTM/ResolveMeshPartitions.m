function [nxParts, nyParts, mode] = ResolveMeshPartitions(lengthXAxis, lengthYAxis, ...
    numPartitionsMicroscale, meshNumPartitionsX, meshNumPartitionsY, meshTargetElementSizeCm)
% ResolveMeshPartitions - Convert mesh-density controls to X/Y partitions.

if nargin < 4
    meshNumPartitionsX = [];
end
if nargin < 5
    meshNumPartitionsY = [];
end
if nargin < 6
    meshTargetElementSizeCm = [];
end

if ~isempty(meshTargetElementSizeCm)
    h = max(eps, meshTargetElementSizeCm(1));
    nxParts = max(1, ceil(lengthXAxis / h));
    nyParts = max(1, ceil(lengthYAxis / h));
    mode = "target_element_size";
elseif ~isempty(meshNumPartitionsX) && ~isempty(meshNumPartitionsY)
    nxParts = max(1, round(meshNumPartitionsX(1)));
    nyParts = max(1, round(meshNumPartitionsY(1)));
    mode = "explicit_xy";
elseif ~isempty(meshNumPartitionsX)
    nxParts = max(1, round(meshNumPartitionsX(1)));
    nyParts = max(1, round(nxParts * lengthYAxis / lengthXAxis));
    mode = "explicit_x";
elseif ~isempty(meshNumPartitionsY)
    nyParts = max(1, round(meshNumPartitionsY(1)));
    nxParts = max(1, round(nyParts * lengthXAxis / lengthYAxis));
    mode = "explicit_y";
else
    aspect = max(1, round(lengthXAxis / lengthYAxis));
    nxParts = max(1, round(numPartitionsMicroscale(1))) * aspect;
    nyParts = max(1, round(numPartitionsMicroscale(1)));
    mode = "legacy_numPartitionsMicroscale";
end
end

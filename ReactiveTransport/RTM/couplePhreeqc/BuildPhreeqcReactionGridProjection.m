function projection = BuildPhreeqcReactionGridProjection(fineGrid, reactionGrid)
% Assign each fine-grid triangle barycenter to one PHREEQC reaction triangle.

fineBary = fineGrid.baryT;
reactionBary = reactionGrid.baryT;
reactionVertices = reactionGrid.coordV;
reactionTriangles = reactionGrid.V0T;

numFine = size(fineBary, 1);
numReaction = size(reactionBary, 1);
fineToReaction = zeros(numFine, 1);

for iFine = 1:numFine
    point = fineBary(iFine, :);
    assigned = 0;
    for iReaction = 1:numReaction
        tri = reactionVertices(reactionTriangles(iReaction, :), :);
        if pointInTriangle(point, tri)
            assigned = iReaction;
            break;
        end
    end
    if assigned == 0
        [~, assigned] = min(sum((reactionBary - point).^2, 2));
    end
    fineToReaction(iFine) = assigned;
end

projection = struct();
projection.fineToReactionCell = fineToReaction;
projection.numFineCells = numFine;
projection.numReactionCells = numReaction;
end

function tf = pointInTriangle(point, tri)
tol = 1e-12;
v0 = tri(3, :) - tri(1, :);
v1 = tri(2, :) - tri(1, :);
v2 = point - tri(1, :);
den = v0(1) * v1(2) - v1(1) * v0(2);
if abs(den) < eps
    tf = false;
    return;
end
a = (v2(1) * v1(2) - v1(1) * v2(2)) / den;
b = (v0(1) * v2(2) - v2(1) * v0(2)) / den;
tf = a >= -tol && b >= -tol && (a + b) <= 1 + tol;
end

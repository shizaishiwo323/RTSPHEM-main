function projection = BuildPhreeqcReactionGridProjection(fineGrid, reactionGrid)
% Assign each fine-grid triangle barycenter to one PHREEQC reaction triangle.

fineBary = fineGrid.baryT;
reactionBary = reactionGrid.baryT;
reactionVertices = reactionGrid.coordV;
reactionTriangles = reactionGrid.V0T;

numFine = size(fineBary, 1);
numReaction = size(reactionBary, 1);
reactionTriangulation = triangulation(reactionTriangles, reactionVertices);
fineToReaction = pointLocation(reactionTriangulation, fineBary);
outside = isnan(fineToReaction);
if any(outside)
    fineToReaction(outside) = dsearchn(reactionBary, fineBary(outside, :));
end

projection = struct();
projection.fineToReactionCell = fineToReaction;
projection.numFineCells = numFine;
projection.numReactionCells = numReaction;
end

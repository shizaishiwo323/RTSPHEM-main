function tests = test_BuildPhreeqcReactionGridProjection
tests = functiontests(localfunctions);
end

function setupOnce(~)
reactiveRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(reactiveRoot));
end

function testMapsFineBarycentersToContainingReactionTriangles(testCase)
fine = domainRectangle(0, 0.06, 0, 0.04, 0.0025);
reaction = domainRectangle(0, 0.06, 0, 0.04, 0.01);

projection = BuildPhreeqcReactionGridProjection(fine, reaction);
tri = triangulation(reaction.V0T, reaction.coordV);
expected = pointLocation(tri, fine.baryT);

verifyEqual(testCase, projection.fineToReactionCell, expected);
verifyEqual(testCase, projection.numFineCells, fine.numT);
verifyEqual(testCase, projection.numReactionCells, reaction.numT);
end

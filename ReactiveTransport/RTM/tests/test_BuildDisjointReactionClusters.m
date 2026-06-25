function tests = test_BuildDisjointReactionClusters
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
rtmDir = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.rtmDir = rtmDir;
addpath(rtmDir);
end

function teardownOnce(~)
% Keep shared MATLAB paths available when directory suites run.
end

function testNonOverlappingCandidatesRemainSeparate(testCase)
sourceCells = [2; 5];
candidateMembers = {[1 2 3], [5 6]};

result = rtm.chemistry.BuildDisjointReactionClusters(sourceCells, candidateMembers, 6);

verifyEqual(testCase, numel(result.clusters), 2);
verifyEqual(testCase, result.clusters(1).source_cells, 2);
verifyEqual(testCase, result.clusters(1).member_cells, [1; 2; 3]);
verifyEqual(testCase, result.clusters(2).source_cells, 5);
verifyEqual(testCase, result.clusters(2).member_cells, [5; 6]);
verifyEqual(testCase, result.membership_count, [1; 1; 1; 0; 1; 1]);
verifyEqual(testCase, max(result.membership_count), 1);
end

function testOverlappingCandidatesAreMerged(testCase)
sourceCells = [2; 4];
candidateMembers = {[1 2 3], [3 4 5]};

result = rtm.chemistry.BuildDisjointReactionClusters(sourceCells, candidateMembers, 5);

verifyEqual(testCase, numel(result.clusters), 1);
verifyEqual(testCase, result.clusters(1).source_cells, [2; 4]);
verifyEqual(testCase, result.clusters(1).member_cells, [1; 2; 3; 4; 5]);
verifyEqual(testCase, result.cluster_id_by_cell, ones(5, 1));
verifyEqual(testCase, max(result.membership_count), 1);
end

function testClusterConstructionIsOrderIndependent(testCase)
forward = rtm.chemistry.BuildDisjointReactionClusters( ...
    [2; 4; 8], {[1 2 3], [3 4 5], [8 9]}, 9);
reverse = rtm.chemistry.BuildDisjointReactionClusters( ...
    [8; 4; 2], {[8 9], [3 4 5], [1 2 3]}, 9);

verifyEqual(testCase, canonicalClusterSets(forward), canonicalClusterSets(reverse));
verifyEqual(testCase, max(forward.membership_count), 1);
verifyEqual(testCase, max(reverse.membership_count), 1);
end

function testInvalidCellIndexIsRejected(testCase)
verifyError(testCase, @() rtm.chemistry.BuildDisjointReactionClusters( ...
    [2], {[1 2 7]}, 6), 'RTSPHEM:Chemistry:InvalidClusterCell');
end

function sets = canonicalClusterSets(result)
sets = strings(numel(result.clusters), 1);
for iCluster = 1:numel(result.clusters)
    sourceText = sprintf('%d,', result.clusters(iCluster).source_cells);
    memberText = sprintf('%d,', result.clusters(iCluster).member_cells);
    sets(iCluster) = "s:" + sourceText + "|m:" + memberText;
end
sets = sort(sets);
end

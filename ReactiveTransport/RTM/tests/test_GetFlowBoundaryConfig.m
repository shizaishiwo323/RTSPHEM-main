function tests = test_GetFlowBoundaryConfig
tests = functiontests(localfunctions);
end

function testBottomToTopBoundaryMapping(testCase)
cfg = GetFlowBoundaryConfig('bottom_to_top', 0.15, 0.266, 10, 20);

verifyEqual(testCase, cfg.inletId, 1);
verifyEqual(testCase, cfg.outletId, 3);
verifyEqual(testCase, cfg.wallIds, [4, 2]);
verifyEqual(testCase, cfg.velocityVector(2.5), [0; 2.5]);
verifyEqual(testCase, cfg.flowLength, 0.266);
verifyEqual(testCase, cfg.crossSectionLength, 0.15);
verifyEqual(testCase, cfg.segmentEdges, linspace(0, 0.266, 6));
verifyEqual(testCase, cfg.axisIndex, 2);
verifyEqual(testCase, cfg.transverseAxisIndex, 1);
end

function testLeftToRightBoundaryMappingIsDefault(testCase)
cfg = GetFlowBoundaryConfig('', 0.15, 0.266, 10, 20);

verifyEqual(testCase, cfg.inletId, 4);
verifyEqual(testCase, cfg.outletId, 2);
verifyEqual(testCase, cfg.wallIds, [1, 3]);
verifyEqual(testCase, cfg.velocityVector(2.5), [2.5; 0]);
verifyEqual(testCase, cfg.flowLength, 0.15);
verifyEqual(testCase, cfg.crossSectionLength, 0.266);
verifyEqual(testCase, cfg.segmentEdges, linspace(0, 0.15, 6));
verifyEqual(testCase, cfg.axisIndex, 1);
verifyEqual(testCase, cfg.transverseAxisIndex, 2);
end

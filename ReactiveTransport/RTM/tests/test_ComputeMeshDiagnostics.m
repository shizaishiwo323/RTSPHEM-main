function tests = test_ComputeMeshDiagnostics
tests = functiontests(localfunctions);
end

function testReportsBasicTriangularMeshStats(testCase)
coord = [0 0; 1 0; 1 1; 0 1];
triangles = [1 2 3; 1 3 4];
levelSet = [-1; -1; 1; -1];

stats = ComputeMeshDiagnostics(coord, triangles, levelSet);

verifyEqual(testCase, stats.num_nodes, 4);
verifyEqual(testCase, stats.num_triangles, 2);
verifyEqual(testCase, stats.num_edges, 5);
verifyEqual(testCase, stats.triangle_area_min, 0.5, 'AbsTol', 1e-12);
verifyEqual(testCase, stats.triangle_area_max, 0.5, 'AbsTol', 1e-12);
verifyEqual(testCase, stats.num_interface_triangles, 2);
verifyEqual(testCase, stats.num_pore_triangles, 0);
verifyEqual(testCase, stats.num_solid_triangles, 0);
verifyGreaterThan(testCase, stats.triangle_quality_min, 0);
verifyLessThanOrEqual(testCase, stats.triangle_quality_max, 1);
end

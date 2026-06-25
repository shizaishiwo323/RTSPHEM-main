function tests = test_FlowFaceMassConservation
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

function testValidatesBalancedTwoCellFaceFluxes(testCase)
flow = twoCellBalancedFlow();

diagnostics = rtm.flow.ValidateFaceFluxDivergence(flow, 2, ...
    struct('absoluteTolerance_cm3_s', 1e-12));

verifyTrue(testCase, diagnostics.accepted);
verifyEqual(testCase, diagnostics.cell_divergence_cm3_s, zeros(2, 1), ...
    'AbsTol', 1e-18);
verifyEqual(testCase, diagnostics.inlet_flow_cm3_s, 1, 'AbsTol', 1e-18);
verifyEqual(testCase, diagnostics.outlet_flow_cm3_s, 1, 'AbsTol', 1e-18);
verifyEqual(testCase, diagnostics.global_residual_cm3_s, 0, 'AbsTol', 1e-18);
verifyEqual(testCase, diagnostics.inlet_outlet_relative_residual, 0, ...
    'AbsTol', 1e-18);
end

function testRejectsMissingOutflowBoundary(testCase)
flow = twoCellBalancedFlow();
flow.boundary_face_cells = 1;
flow.boundary_face_area_cm2 = 1;
flow.boundary_face_velocity_cm_s = 1;
flow.boundary_type = "dirichlet";

diagnostics = rtm.flow.ValidateFaceFluxDivergence(flow, 2, ...
    struct('absoluteTolerance_cm3_s', 1e-12));

verifyFalse(testCase, diagnostics.accepted);
verifyEqual(testCase, diagnostics.cell_divergence_cm3_s, [0; 1], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, diagnostics.global_residual_cm3_s, 1, 'AbsTol', 1e-18);
verifyTrue(testCase, any(diagnostics.failure_reasons == ...
    "cell divergence exceeds tolerance"));
end

function testRejectsInletOutletRelativeImbalance(testCase)
flow = twoCellBalancedFlow();
flow.boundary_face_velocity_cm_s = [1; 0.99];

diagnostics = rtm.flow.ValidateFaceFluxDivergence(flow, 2, ...
    struct('absoluteTolerance_cm3_s', 1e6, ...
    'relativeTolerance', 1e-3));

verifyFalse(testCase, diagnostics.accepted);
verifyEqual(testCase, diagnostics.inlet_flow_cm3_s, 1, 'AbsTol', 1e-18);
verifyEqual(testCase, diagnostics.outlet_flow_cm3_s, 0.99, 'AbsTol', 1e-18);
verifyEqual(testCase, diagnostics.inlet_outlet_relative_residual, 0.01, ...
    'RelTol', 1e-12, 'AbsTol', 1e-18);
verifyTrue(testCase, any(diagnostics.failure_reasons == ...
    "inlet/outlet flow imbalance exceeds tolerance"));
end

function testMapsRt0FaceFluxesToCutCellTransportOptions(testCase)
hyphm = struct();
hyphm.internal_face_cells = [1 2; 2 3];
hyphm.internal_face_flux_cm3_s = [0.25; -0.5];
hyphm.internal_face_area_cm2 = [2; 4];
hyphm.internal_face_distance_cm = [0.1; 0.2];
hyphm.boundary_face_cells = [1; 3];
hyphm.boundary_face_flux_cm3_s = [0.25; 0.5];
hyphm.boundary_face_area_cm2 = [1; 2];
hyphm.boundary_type = ["dirichlet"; "outflow"];

flow = rtm.flow.MapRt0FluxToCutCellFaces(hyphm);

verifyEqual(testCase, flow.internal_face_cells, hyphm.internal_face_cells);
verifyEqual(testCase, flow.internal_face_area_cm2, [2; 4], 'AbsTol', 1e-18);
verifyEqual(testCase, flow.internal_face_distance_cm, [0.1; 0.2], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, flow.internal_face_velocity_cm_s, [0.125; -0.125], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, flow.boundary_face_velocity_cm_s, [0.25; 0.25], ...
    'AbsTol', 1e-18);
verifyEqual(testCase, flow.boundary_type, ["dirichlet"; "outflow"]);
end

function testBuildFaceFluxesFromHyphmKeepsTransportContract(testCase)
hyphm = struct();
hyphm.internal_face_cells = [1 2];
hyphm.internal_face_flux_cm3_s = 1;
hyphm.internal_face_area_cm2 = 2;
hyphm.internal_face_distance_cm = 0.5;
hyphm.boundary_face_cells = [1; 2];
hyphm.boundary_face_flux_cm3_s = [1; 1];
hyphm.boundary_face_area_cm2 = [1; 1];
hyphm.boundary_type = ["dirichlet"; "outflow"];

flow = rtm.flow.BuildFaceFluxesFromHyphm(hyphm, 2, ...
    struct('validateDivergence', true, 'absoluteTolerance_cm3_s', 1e-12));

verifyTrue(testCase, flow.diagnostics.accepted);
verifyEqual(testCase, flow.internal_face_velocity_cm_s, 0.5, ...
    'AbsTol', 1e-18);
verifyTrue(testCase, isfield(flow, 'boundary_face_cells'));
end

function testConvertsFlowToCutCellTransportOptions(testCase)
flow = twoCellBalancedFlow();
options = struct('diffusion_coefficient_cm2_s', [1e-5; 1e-5]);

transportOptions = rtm.flow.ApplyFaceFluxesToTransportOptions(options, flow);

verifyEqual(testCase, transportOptions.internal_face_cells, ...
    flow.internal_face_cells);
verifyEqual(testCase, transportOptions.internal_face_velocity_cm_s, ...
    flow.internal_face_velocity_cm_s);
verifyEqual(testCase, transportOptions.boundary_face_cells, ...
    flow.boundary_face_cells);
verifyEqual(testCase, transportOptions.boundary_type, flow.boundary_type);
verifyEqual(testCase, transportOptions.diffusion_coefficient_cm2_s, [1e-5; 1e-5]);
end

function testSolveHyphmStokesLevelCallsConfiguredSolverAndMapsFlux(testCase)
captured = struct();
geometry = struct('water_volume_cm3', [1; 1]);
options = struct();
options.solverFunction = @mockSolver;
options.validateDivergence = true;
options.absoluteTolerance_cm3_s = 1e-12;

flow = rtm.flow.SolveHyphmStokesLevel(geometry, options);

verifyEqual(testCase, captured.geometry.water_volume_cm3, [1; 1]);
verifyTrue(testCase, captured.options.validateDivergence);
verifyEqual(testCase, flow.internal_face_velocity_cm_s, 1, 'AbsTol', 1e-18);
verifyTrue(testCase, flow.diagnostics.accepted);

    function hyphm = mockSolver(inputGeometry, inputOptions)
        captured.geometry = inputGeometry;
        captured.options = inputOptions;
        hyphm = twoCellBalancedHyphmFlux();
    end
end

function testSolveHyphmStokesLevelErrorsWithoutSolver(testCase)
verifyError(testCase, ...
    @() rtm.flow.SolveHyphmStokesLevel(struct(), struct()), ...
    'RTSPHEM:Flow:MissingHyphmSolver');
end

function flow = twoCellBalancedFlow()
flow = struct();
flow.internal_face_cells = [1 2];
flow.internal_face_area_cm2 = 1;
flow.internal_face_distance_cm = 1;
flow.internal_face_velocity_cm_s = 1;
flow.boundary_face_cells = [1; 2];
flow.boundary_face_area_cm2 = [1; 1];
flow.boundary_face_velocity_cm_s = [1; 1];
flow.boundary_type = ["dirichlet"; "outflow"];
end

function hyphm = twoCellBalancedHyphmFlux()
hyphm = struct();
hyphm.internal_face_cells = [1 2];
hyphm.internal_face_flux_cm3_s = 1;
hyphm.internal_face_area_cm2 = 1;
hyphm.internal_face_distance_cm = 1;
hyphm.boundary_face_cells = [1; 2];
hyphm.boundary_face_flux_cm3_s = [1; 1];
hyphm.boundary_face_area_cm2 = [1; 1];
hyphm.boundary_type = ["dirichlet"; "outflow"];
end

function tests = test_NMRSimulationConfig
tests = functiontests(localfunctions);
end

function testDefaultConfigHasUnifiedMethod(testCase)
cfg = NMRSimulationConfig();
verifyTrue(testCase, isfield(cfg, 'nmr_method'));
verifyTrue(testCase, isfield(cfg, 'nmrSurrogateModelPath'));
verifyTrue(testCase, isfield(cfg, 'mph_file'));
verifyTrue(testCase, isfield(cfg, 'python_exe'));
end

function testUnifiedMethodSelection(testCase)
opts = ResolveNMRSimulationOptions(struct('nmr_method', 'comsol'));
verifyTrue(testCase, opts.enableNMRSimulation);
verifyFalse(testCase, opts.enableNMRSurrogate);
verifyFalse(testCase, opts.enablePNGSimulation);

opts = ResolveNMRSimulationOptions(struct('nmr_method', 'surrogate'));
verifyFalse(testCase, opts.enableNMRSimulation);
verifyTrue(testCase, opts.enableNMRSurrogate);
verifyFalse(testCase, opts.enablePNGSimulation);

opts = ResolveNMRSimulationOptions(struct('nmr_method', 'png_mesh'));
verifyFalse(testCase, opts.enableNMRSimulation);
verifyFalse(testCase, opts.enableNMRSurrogate);
verifyTrue(testCase, opts.enablePNGSimulation);
verifyEqual(testCase, opts.pngNMRMethod, 'png_mesh');

opts = ResolveNMRSimulationOptions(struct('nmr_method', 'none'));
verifyFalse(testCase, opts.enableNMRSimulation);
verifyFalse(testCase, opts.enableNMRSurrogate);
verifyFalse(testCase, opts.enablePNGSimulation);
end

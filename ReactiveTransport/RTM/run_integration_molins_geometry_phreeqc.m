function run = run_integration_molins_geometry_phreeqc(mode, options)
%RUN_INTEGRATION_MOLINS_GEOMETRY_PHREEQC Run the Molins PHREEQC integration suite.
%
% By default this entry uses the real PHREEQC path configured by the
% benchmark case factory. Tests and dry runs may pass
% options.phreeqcRunBatchFunction to inject a mock PHREEQC batch runner.

if nargin < 1 || isempty(mode)
    mode = 'integration';
end
if nargin < 2 || isempty(options)
    options = struct();
end
options.kinds = "integration_phreeqc";
run = run_molins_convergence_suites(mode, options);
end

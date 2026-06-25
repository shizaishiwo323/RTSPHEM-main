function run = run_benchmark_molins_partII_strict(mode, options)
%RUN_BENCHMARK_MOLINS_PARTII_STRICT Run the strict Molins Part II suite.

if nargin < 1 || isempty(mode)
    mode = 'diagnostic';
end
if nargin < 2 || isempty(options)
    options = struct();
end
options.kinds = "partII_strict";
run = run_molins_convergence_suites(mode, options);
end

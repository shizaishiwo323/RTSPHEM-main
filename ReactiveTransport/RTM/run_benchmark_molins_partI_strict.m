function run = run_benchmark_molins_partI_strict(mode, options)
%RUN_BENCHMARK_MOLINS_PARTI_STRICT Run the strict Molins Part I suite.

if nargin < 1 || isempty(mode)
    mode = 'diagnostic';
end
if nargin < 2 || isempty(options)
    options = struct();
end
options.kinds = "partI_strict";
run = run_molins_convergence_suites(mode, options);
end

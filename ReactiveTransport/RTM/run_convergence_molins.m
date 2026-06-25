function run = run_convergence_molins(mode, options)
%RUN_CONVERGENCE_MOLINS Compatibility entry for Molins convergence suites.

if nargin < 1
    run = run_molins_convergence_suites();
elseif nargin < 2
    run = run_molins_convergence_suites(mode);
else
    run = run_molins_convergence_suites(mode, options);
end
end

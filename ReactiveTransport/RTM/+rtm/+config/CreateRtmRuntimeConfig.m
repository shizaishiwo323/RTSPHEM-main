function runtimeConfig = CreateRtmRuntimeConfig(config, manifestOverrides)
%CREATERTMRUNTIMECONFIG Validate RTM config and attach runtime provenance.

if nargin < 1 || isempty(config)
    config = struct();
end
if nargin < 2 || isempty(manifestOverrides)
    manifestOverrides = struct();
end

validatedConfig = rtm.config.ValidateReactiveTransportConfig(config);
rtmDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
phreeqcRuntime = rtm.phreeqc.ResolveRuntime(rtmDir, validatedConfig);
manifest = rtm.diagnostics.CreateRuntimeManifest(validatedConfig, ...
    manifestOverrides);

runtimeConfig = struct();
runtimeConfig.config = validatedConfig;
runtimeConfig.phreeqc_runtime = phreeqcRuntime;
runtimeConfig.runtime_manifest = manifest;
runtimeConfig.solver_architecture = manifest.solver_architecture;
runtimeConfig.chemistry_mode = manifest.chemistry_mode;
runtimeConfig.transport_backend = manifest.transport_backend;
runtimeConfig.operator_order = manifest.operator_order;
runtimeConfig.units = manifest.units;
end

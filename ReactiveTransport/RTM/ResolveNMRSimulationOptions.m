function options = ResolveNMRSimulationOptions(config)
% ResolveNMRSimulationOptions - Convert unified NMR method config to switches.
%
% Supported methods:
%   none, off, disabled
%   comsol
%   surrogate, ml, machine_learning, nmr_agent
%   png_pixel_cpu, png_pixel_gpu, png_mesh

if nargin < 1 || isempty(config)
    config = NMRSimulationConfig();
end

method = char(cfgget(config, 'nmr_method', cfgget(config, 'nmrMethod', '')));
if isempty(strtrim(method))
    method = inferLegacyMethod(config);
end
method = lower(strtrim(method));

options = struct();
options.nmr_method = method;
options.enableNMRSimulation = false;
options.enableNMRSurrogate = false;
options.enablePNGSimulation = false;
options.pngNMRMethod = '';
options.config = config;

switch method
    case {'none', 'off', 'disabled', 'false', 'no'}
        return;
    case {'comsol', 'comsol_nmr'}
        options.enableNMRSimulation = true;
    case {'surrogate', 'ml', 'machine_learning', 'nmr_agent', 'nmr-agent'}
        options.enableNMRSurrogate = true;
    case {'png_pixel_cpu', 'png_pixel_gpu', 'png_mesh'}
        options.enablePNGSimulation = true;
        options.pngNMRMethod = method;
        options.config.method = method;
    otherwise
        error('NMRConfig:UnsupportedMethod', ...
            'Unsupported nmr_method "%s". Use none, comsol, surrogate, png_pixel_cpu, png_pixel_gpu, or png_mesh.', ...
            method);
end
end

function method = inferLegacyMethod(config)
if logical(cfgget(config, 'enableNMRSimulation', false))
    method = 'comsol';
elseif logical(cfgget(config, 'enableNMRSurrogate', false))
    method = 'surrogate';
elseif logical(cfgget(config, 'enablePNGSimulation', false))
    method = char(cfgget(config, 'pngNMRMethod', cfgget(config, 'method', 'png_mesh')));
else
    method = 'none';
end
end

function value = cfgget(config, fieldName, defaultValue)
if isstruct(config) && isfield(config, fieldName) && ~isempty(config.(fieldName))
    value = config.(fieldName);
else
    value = defaultValue;
end
end

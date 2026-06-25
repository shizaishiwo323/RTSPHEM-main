function [stableDt, diagnostics] = precip_ComputeTransportStableDt(spec, options)
% precip_ComputeTransportStableDt - CFL-based finite-volume transport step.
%
% Inputs:
%   spec    - benchmark spec with dx, dy, velocity, and diffusivity.
%   options - advectiveCfl and diffusiveCfl controls.
%
% Output:
%   stableDt    - minimum positive explicit transport stability step, seconds.
%   diagnostics - individual limits and active limiter label.

if nargin < 2
    options = struct();
end

advectiveCfl = getFieldOrDefault(options, 'advectiveCfl', 0.5);
diffusiveCfl = getFieldOrDefault(options, 'diffusiveCfl', 0.25);
if advectiveCfl <= 0 || diffusiveCfl <= 0
    error('RTSPHEM:Precipitate:InvalidTransportCfl', ...
        'Transport CFL controls must be positive.');
end

u = abs(getFieldOrDefault(spec, 'darcyVelocity_cm_s', 0));
flowField = getFieldOrDefault(options, 'flowField', []);
if ~isempty(flowField)
    [u, maxV] = maxVelocityFromFlowField(flowField, spec);
else
    maxV = 0;
end
d = maxDiffusivity(spec, options);
if u > 0 || maxV > 0
    xLimit = Inf;
    yLimit = Inf;
    if u > 0
        xLimit = spec.dx_cm / u;
    end
    if maxV > 0
        yLimit = spec.dy_cm / maxV;
    end
    advectiveLimit = advectiveCfl * min(xLimit, yLimit);
else
    advectiveLimit = Inf;
end
if d > 0
    diffusiveLimit = diffusiveCfl / (2 * d * ...
        (1 / spec.dx_cm^2 + 1 / spec.dy_cm^2));
else
    diffusiveLimit = Inf;
end

stableDt = min(advectiveLimit, diffusiveLimit);
if ~isfinite(stableDt)
    stableDt = Inf;
    limiter = "none";
elseif advectiveLimit <= diffusiveLimit
    limiter = "advective_cfl";
else
    limiter = "diffusive_cfl";
end

diagnostics = struct();
diagnostics.advectiveLimit_s = advectiveLimit;
diagnostics.diffusiveLimit_s = diffusiveLimit;
diagnostics.stableDt_s = stableDt;
diagnostics.limiter = limiter;
diagnostics.advectiveCfl = advectiveCfl;
diagnostics.diffusiveCfl = diffusiveCfl;
diagnostics.maxVelocity_cm_s = max(u, maxV);
diagnostics.maxDiffusivity_cm2_s = d;
end

function [maxU, maxV] = maxVelocityFromFlowField(flowField, spec)
if ~isfield(flowField, 'velocityX_cm_s') || isempty(flowField.velocityX_cm_s)
    error('RTSPHEM:Precipitate:MissingVelocityField', ...
        'flowField.velocityX_cm_s is required for flow-field CFL.');
end
u = flowField.velocityX_cm_s;
if ~isequal(size(u), [spec.numY, spec.numX])
    error('RTSPHEM:Precipitate:InvalidVelocityFieldSize', ...
        'flowField.velocityX_cm_s must have size [numY, numX].');
end
if isfield(flowField, 'velocityY_cm_s') && ~isempty(flowField.velocityY_cm_s)
    v = flowField.velocityY_cm_s;
    if ~isequal(size(v), [spec.numY, spec.numX])
        error('RTSPHEM:Precipitate:InvalidVelocityFieldSize', ...
            'flowField.velocityY_cm_s must have size [numY, numX].');
    end
else
    v = zeros(size(u));
end
if any(~isfinite(u(:))) || any(~isfinite(v(:)))
    error('RTSPHEM:Precipitate:InvalidVelocityFieldValue', ...
        'Flow-field velocities must be finite.');
end
maxU = max(abs(u(:)));
maxV = max(abs(v(:)));
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function d = maxDiffusivity(spec, options)
if isfield(options, 'effectiveDiffusivity_cm2_s') && ...
        ~isempty(options.effectiveDiffusivity_cm2_s)
    dField = options.effectiveDiffusivity_cm2_s;
elseif isfield(options, 'effectiveDiffusivity') && ...
        ~isempty(options.effectiveDiffusivity)
    dField = options.effectiveDiffusivity;
else
    d = max(getFieldOrDefault(spec, 'diffusionCoefficient_cm2_s', 0), 0);
    return;
end
if ~isequal(size(dField), [spec.numY, spec.numX])
    error('RTSPHEM:Precipitate:InvalidDiffusivityFieldSize', ...
        'effectiveDiffusivity_cm2_s must have size [numY, numX].');
end
if any(~isfinite(dField(:))) || any(dField(:) < 0)
    error('RTSPHEM:Precipitate:InvalidDiffusivityFieldValue', ...
        'effectiveDiffusivity_cm2_s must contain finite nonnegative values.');
end
d = max(dField(:));
end

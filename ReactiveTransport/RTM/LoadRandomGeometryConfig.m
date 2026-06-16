function [loadedData, params] = LoadRandomGeometryConfig(configPath, defaults)
% LoadRandomGeometryConfig - Load and normalize a saved random PNM geometry.
%
% Inputs:
%   configPath - path to random_geometry_config.mat
%   defaults   - current caller parameters; loaded geometry fields override them
%
% Outputs:
%   loadedData - raw MAT variables loaded from configPath
%   params     - geometry parameters to use for the run

if nargin < 2 || isempty(defaults)
    defaults = struct();
end

params = defaults;
configPath = char(configPath);
if isempty(strtrim(configPath))
    error('RTM:RandomGeometryConfig:EmptyPath', '随机几何配置路径为空。');
end
if exist(configPath, 'file') ~= 2
    error('RTM:RandomGeometryConfig:MissingFile', ...
        '未找到随机几何配置文件: %s', configPath);
end
if isGitLfsPointerFile(configPath)
    error('RTM:RandomGeometryConfig:GitLfsPointer', ...
        ['随机几何配置文件是 Git LFS pointer，不是真正的 MAT 文件: %s\n', ...
         '请先运行 git lfs pull 或从原始数据位置取回实际 random_geometry_config.mat。'], ...
        configPath);
end

loadedData = load(configPath);
validateLoadedGeometry(loadedData, configPath);

params.circleRadius = getFirstNumericField(loadedData, ...
    {'circleRadius', 'particleRadius_cm'}, getDefault(params, 'circleRadius', []));
if isempty(params.circleRadius) && isfield(loadedData, 'circleRadii')
    params.circleRadius = mean(loadedData.circleRadii(:));
end

params.circleSpacing = getFirstNumericField(loadedData, ...
    {'circleSpacing'}, getDefault(params, 'circleSpacing', []));
params.targetAvgSpacing = getFirstNumericField(loadedData, ...
    {'targetAvgSpacing', 'finalAvgSpacing'}, getDefault(params, 'targetAvgSpacing', []));
params.minThroatRandom = getFirstNumericField(loadedData, ...
    {'minThroatRandom'}, getDefault(params, 'minThroatRandom', []));

params.lengthXAxis = loadedData.lengthXAxis;
params.lengthYAxis = loadedData.lengthYAxis;
params.targetLengthYAxis = loadedData.lengthYAxis;
params.targetAspectRatio = loadedData.lengthXAxis / loadedData.lengthYAxis;

if isfield(loadedData, 'circleRadii') && ~isempty(loadedData.circleRadii)
    radii = loadedData.circleRadii(:);
    loadedData.circleRadii = radii;
    params.useRandomParticleRadii = any(abs(radii - radii(1)) > 1e-12);
    params.randomParticleRadiusMin = min(radii);
    params.randomParticleRadiusMax = max(radii);
else
    params.useRandomParticleRadii = getDefault(params, 'useRandomParticleRadii', false);
    params.randomParticleRadiusMin = getDefault(params, 'randomParticleRadiusMin', params.circleRadius);
    params.randomParticleRadiusMax = getDefault(params, 'randomParticleRadiusMax', params.circleRadius);
end

if isfield(loadedData, 'targetInitialPorosity') && ~isempty(loadedData.targetInitialPorosity)
    params.targetInitialPorosity = loadedData.targetInitialPorosity;
end
end

function validateLoadedGeometry(loadedData, configPath)
requiredFields = {'circleCenters', 'lengthXAxis', 'lengthYAxis'};
for iField = 1:numel(requiredFields)
    fieldName = requiredFields{iField};
    if ~isfield(loadedData, fieldName) || isempty(loadedData.(fieldName))
        error('RTM:RandomGeometryConfig:MissingField', ...
            '随机几何配置缺少字段 %s: %s', fieldName, configPath);
    end
end
if ~isnumeric(loadedData.circleCenters) || size(loadedData.circleCenters, 2) ~= 2
    error('RTM:RandomGeometryConfig:InvalidCenters', ...
        'circleCenters 必须是 N×2 数值数组: %s', configPath);
end
if ~isfiniteScalar(loadedData.lengthXAxis) || loadedData.lengthXAxis <= 0 || ...
        ~isfiniteScalar(loadedData.lengthYAxis) || loadedData.lengthYAxis <= 0
    error('RTM:RandomGeometryConfig:InvalidSize', ...
        'lengthXAxis/lengthYAxis 必须为正数: %s', configPath);
end
if isfield(loadedData, 'circleRadii') && ~isempty(loadedData.circleRadii)
    if numel(loadedData.circleRadii) ~= size(loadedData.circleCenters, 1)
        error('RTM:RandomGeometryConfig:InvalidRadii', ...
            'circleRadii 数量必须与 circleCenters 行数一致: %s', configPath);
    end
elseif ~isfield(loadedData, 'circleRadius') || isempty(loadedData.circleRadius)
    error('RTM:RandomGeometryConfig:MissingRadius', ...
        '随机几何配置需要 circleRadius 或 circleRadii: %s', configPath);
end
end

function tf = isGitLfsPointerFile(configPath)
fid = fopen(configPath, 'r');
if fid < 0
    tf = false;
    return;
end
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, 256, '*char')';
tf = startsWith(bytes, 'version https://git-lfs.github.com/spec/v1');
end

function value = getFirstNumericField(s, names, fallback)
value = fallback;
for iName = 1:numel(names)
    name = names{iName};
    if isfield(s, name) && ~isempty(s.(name)) && isnumeric(s.(name))
        value = s.(name);
        value = value(1);
        return;
    end
end
end

function value = getDefault(s, name, fallback)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = fallback;
end
end

function tf = isfiniteScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value);
end

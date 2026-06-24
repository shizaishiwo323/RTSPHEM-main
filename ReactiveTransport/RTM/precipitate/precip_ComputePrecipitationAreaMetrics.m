function metrics = precip_ComputePrecipitationAreaMetrics(grid, levelSetNow, levelSetInitial, config)
% precip_ComputePrecipitationAreaMetrics - Compute solid-area metrics for precipitation benchmarks.
%
% Window metrics use approximate x-coordinate windows for the local
% Zhang/Yoon-style cylindrical-post geometry. Strict paper-image
% reproduction should replace these windows with digitized geometry masks.

if nargin < 4 || isempty(config)
    config = struct();
end

levelSetNow = levelSetNow(:);
levelSetInitial = levelSetInitial(:);
validateGridAndLevels(grid, levelSetNow, levelSetInitial);

areaConfig = configureAreaWindows(config);
totalMask = true(numel(grid.areaT), 1);
firstPoreMask = grid.baryT(:, 1) <= areaConfig.firstPoreXMaxCm;
firstThreePoresMask = grid.baryT(:, 1) <= areaConfig.firstThreePoresXMaxCm;

currentTotal = computeSolidAreaForTriangles(grid, levelSetNow, totalMask);
initialTotal = computeSolidAreaForTriangles(grid, levelSetInitial, totalMask);
currentFirst = computeSolidAreaForTriangles(grid, levelSetNow, firstPoreMask);
initialFirst = computeSolidAreaForTriangles(grid, levelSetInitial, firstPoreMask);
currentFirstThree = computeSolidAreaForTriangles(grid, levelSetNow, firstThreePoresMask);
initialFirstThree = computeSolidAreaForTriangles(grid, levelSetInitial, firstThreePoresMask);

metrics = struct();
metrics.totalSolidArea_cm2 = currentTotal;
metrics.totalNetSolidArea_cm2 = currentTotal - initialTotal;
metrics.firstPoreSolidArea_cm2 = currentFirst;
metrics.firstPoreNetSolidArea_cm2 = currentFirst - initialFirst;
metrics.firstThreePoresSolidArea_cm2 = currentFirstThree;
metrics.firstThreePoresNetSolidArea_cm2 = currentFirstThree - initialFirstThree;
metrics.firstPoreXMaxCm = areaConfig.firstPoreXMaxCm;
metrics.firstThreePoresXMaxCm = areaConfig.firstThreePoresXMaxCm;
metrics.windowDefinition = areaConfig.windowDefinition;
end

function validateGridAndLevels(grid, levelSetNow, levelSetInitial)
requiredFields = {'coordV', 'V0T', 'areaT', 'baryT'};
for iField = 1:numel(requiredFields)
    if ~isstruct(grid) && ~isobject(grid)
        error('RTSPHEM:Precipitate:InvalidAreaGrid', ...
            'Area metrics require a grid struct/object.');
    end
    if ~isfield_or_prop(grid, requiredFields{iField})
        error('RTSPHEM:Precipitate:InvalidAreaGrid', ...
            'Area metrics grid is missing %s.', requiredFields{iField});
    end
end
numVertices = size(grid.coordV, 1);
numTriangles = size(grid.V0T, 1);
if ~isnumeric(grid.coordV) || size(grid.coordV, 2) ~= 2 || ...
        isempty(grid.coordV) || any(~isfinite(grid.coordV(:)))
    error('RTSPHEM:Precipitate:InvalidAreaGrid', ...
        'Grid coordV must be a finite N-by-2 coordinate array.');
end
if ~isnumeric(grid.V0T) || size(grid.V0T, 2) ~= 3 || ...
        isempty(grid.V0T) || any(~isfinite(grid.V0T(:))) || ...
        any(abs(grid.V0T(:) - round(grid.V0T(:))) > 0) || ...
        any(grid.V0T(:) < 1) || any(grid.V0T(:) > numVertices)
    error('RTSPHEM:Precipitate:InvalidAreaGrid', ...
        'Grid V0T must be a finite integer N-by-3 array of valid vertex indices.');
end
if ~isnumeric(grid.areaT) || numel(grid.areaT) ~= numTriangles || ...
        any(~isfinite(grid.areaT(:))) || any(grid.areaT(:) < 0)
    error('RTSPHEM:Precipitate:InvalidAreaGrid', ...
        'Grid areaT must be finite, nonnegative, and have one value per triangle.');
end
if ~isnumeric(grid.baryT) || size(grid.baryT, 1) ~= numTriangles || ...
        size(grid.baryT, 2) < 2 || any(~isfinite(grid.baryT(:)))
    error('RTSPHEM:Precipitate:InvalidAreaGrid', ...
        'Grid baryT must be finite and have at least x/y columns for every triangle.');
end
if ~isnumeric(levelSetNow) || ~isreal(levelSetNow) || ...
        ~isnumeric(levelSetInitial) || ~isreal(levelSetInitial) || ...
        numel(levelSetNow) ~= numVertices || numel(levelSetInitial) ~= numVertices || ...
        any(~isfinite(levelSetNow(:))) || any(~isfinite(levelSetInitial(:)))
    error('RTSPHEM:Precipitate:InvalidAreaLevels', ...
        'Level-set arrays must be finite real numeric values, one per grid vertex.');
end
end

function tf = isfield_or_prop(value, name)
tf = (isstruct(value) && isfield(value, name)) || ...
    (isobject(value) && isprop(value, name));
end

function areaConfig = configureAreaWindows(config)
postDiameter = cfgget(config, 'postDiameter', 0.03);
poreBody = cfgget(config, 'poreBody', 0.018);
poreThroat = cfgget(config, 'poreThroat', 0.004);
inletClearance = cfgget(config, 'zhangInletClearanceCm', max(poreBody, poreThroat));
pitchX = postDiameter + poreBody;
defaultFirstPoreXMax = inletClearance + postDiameter + poreBody + poreThroat;
defaultFirstThreeXMax = inletClearance + 3 * pitchX + poreThroat;

areaConfig = struct();
areaConfig.firstPoreXMaxCm = cfgget(config, 'benchmarkFirstPoreXMaxCm', defaultFirstPoreXMax);
areaConfig.firstThreePoresXMaxCm = cfgget(config, ...
    'benchmarkFirstThreePoresXMaxCm', defaultFirstThreeXMax);
areaConfig.windowDefinition = 'approximate_x_windows';
end

function solidArea = computeSolidAreaForTriangles(grid, levels, triangleMask)
solidArea = 0;
triangleIds = find(triangleMask(:));
for iTriangle = 1:numel(triangleIds)
    kT = triangleIds(iTriangle);
    vertexIds = grid.V0T(kT, :);
    points = grid.coordV(vertexIds, :);
    values = levels(vertexIds);
    waterArea = triangleAreaBelowLevelZero(points, values);
    solidArea = solidArea + max(grid.areaT(kT) - waterArea, 0);
end
end

function area = triangleAreaBelowLevelZero(points, values)
inputPoints = points;
inputValues = values(:);
polyPoints = zeros(0, 2);
nInput = size(inputPoints, 1);
for iPoint = 1:nInput
    jPoint = mod(iPoint, nInput) + 1;
    p1 = inputPoints(iPoint, :);
    p2 = inputPoints(jPoint, :);
    v1 = inputValues(iPoint);
    v2 = inputValues(jPoint);
    inside1 = v1 <= 0;
    inside2 = v2 <= 0;
    if inside1 && inside2
        polyPoints(end + 1, :) = p2; %#ok<AGROW>
    elseif inside1 && ~inside2
        polyPoints(end + 1, :) = interpolateLevelCrossing(p1, p2, v1, v2); %#ok<AGROW>
    elseif ~inside1 && inside2
        polyPoints(end + 1, :) = interpolateLevelCrossing(p1, p2, v1, v2); %#ok<AGROW>
        polyPoints(end + 1, :) = p2; %#ok<AGROW>
    end
end

if size(polyPoints, 1) < 3
    area = 0;
else
    x = polyPoints(:, 1);
    y = polyPoints(:, 2);
    area = 0.5 * abs(sum(x .* y([2:end, 1]) - y .* x([2:end, 1])));
end
end

function point = interpolateLevelCrossing(p1, p2, v1, v2)
denominator = v1 - v2;
if abs(denominator) < eps
    lambda = 0.5;
else
    lambda = v1 / denominator;
end
lambda = min(max(lambda, 0), 1);
point = p1 + lambda * (p2 - p1);
end

function value = cfgget(config, fieldName, defaultValue)
if isstruct(config) && isfield(config, fieldName) && ~isempty(config.(fieldName))
    value = config.(fieldName);
else
    value = defaultValue;
end
end

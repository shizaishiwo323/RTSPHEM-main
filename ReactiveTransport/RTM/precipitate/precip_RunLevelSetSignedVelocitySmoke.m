function smoke = precip_RunLevelSetSignedVelocitySmoke(options)
% precip_RunLevelSetSignedVelocitySmoke - Verify signed level-set speed convention.
%
% The precipitation solver uses levelSet > 0 for solid and levelSet < 0 for
% pore. Positive normal speed is the legacy dissolution direction, so a
% signed precipitation speed is negative and must increase solid area.
% Solid area is estimated by node counting for this lightweight smoke test;
% production benchmark area metrics should use triangle/subcell geometry.

if nargin < 1
    options = struct();
end

domainLengthCm = cfgget(options, 'domainLengthCm', 1.0);
numSamples = cfgget(options, 'numSamples', 401);
numPartitions = numSamples - 1;
initialRadiusCm = cfgget(options, 'initialRadiusCm', 0.20);
timeStepSize = cfgget(options, 'timeStepSize', 1.0);
speedMagnitudeCmS = cfgget(options, 'speedMagnitudeCmS', 0.02);

if numSamples < 21
    error('RTSPHEM:Precipitate:InvalidSmokeGrid', ...
        'numSamples must be at least 21 for the level-set smoke test.');
end
if initialRadiusCm <= 0 || initialRadiusCm >= 0.5 * domainLengthCm
    error('RTSPHEM:Precipitate:InvalidSmokeRadius', ...
        'initialRadiusCm must fit inside the square smoke-test domain.');
end
if timeStepSize <= 0 || speedMagnitudeCmS <= 0
    error('RTSPHEM:Precipitate:InvalidSmokeStep', ...
        'timeStepSize and speedMagnitudeCmS must be positive.');
end

grid = FoldedCartesianGrid(2, [0, domainLengthCm, 0, domainLengthCm], ...
    [numPartitions, numPartitions]);
xCoord = grid.coordinates(:, 1);
yCoord = grid.coordinates(:, 2);
center = 0.5 * domainLengthCm;
distanceFromCenter = hypot(xCoord - center, yCoord - center);

initialLevelSet = initialRadiusCm - distanceFromCenter;
precipitationNormalSpeed = -speedMagnitudeCmS;
dissolutionNormalSpeed = speedMagnitudeCmS;

precipitationLevelSet = advanceSignedLevelSet(grid, initialLevelSet, ...
    precipitationNormalSpeed, timeStepSize);
dissolutionLevelSet = advanceSignedLevelSet(grid, initialLevelSet, ...
    dissolutionNormalSpeed, timeStepSize);

cellAreaCm2 = prod(grid.stepSize);
smoke = struct();
smoke.initialSolidArea_cm2 = solidAreaFromLevelSet(initialLevelSet, cellAreaCm2);
smoke.precipitationSolidAreaAfter_cm2 = solidAreaFromLevelSet(precipitationLevelSet, cellAreaCm2);
smoke.dissolutionSolidAreaAfter_cm2 = solidAreaFromLevelSet(dissolutionLevelSet, cellAreaCm2);
smoke.precipitationNormalSpeed_cm_s = precipitationNormalSpeed;
smoke.dissolutionNormalSpeed_cm_s = dissolutionNormalSpeed;
smoke.timeStepSize_s = timeStepSize;
smoke.initialRadius_cm = initialRadiusCm;
smoke.precipitationExpectedRadius_cm = initialRadiusCm - precipitationNormalSpeed * timeStepSize;
smoke.dissolutionExpectedRadius_cm = initialRadiusCm - dissolutionNormalSpeed * timeStepSize;
smoke.levelSetStepper = 'levelSetEquationTimeStep';
end

function levelSetAfter = advanceSignedLevelSet(grid, levelSetBefore, normalSpeed, timeStepSize)
levelSetAfter = levelSetEquationTimeStep(timeStepSize, 0, levelSetBefore, ...
    grid, normalSpeed, 1);
end

function area = solidAreaFromLevelSet(levelSet, cellAreaCm2)
area = nnz(levelSet > 0) * cellAreaCm2;
end

function value = cfgget(config, fieldName, defaultValue)
if isstruct(config) && isfield(config, fieldName) && ~isempty(config.(fieldName))
    value = config.(fieldName);
else
    value = defaultValue;
end
end

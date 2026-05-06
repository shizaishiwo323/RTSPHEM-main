function result = PNM_batch(params, batchResultsDir, expIdx, permeabilityRatioThreshold)
% PNM_batch - Batch adapter for the canonical PNM_beauty3 solver.
%
% This file intentionally does not contain a second copy of the reactive
% transport solver. It converts batch parameters into the config struct used by
% PNM_beauty3, so single-run and batch-run behavior stay aligned.

if nargin < 4 || isempty(permeabilityRatioThreshold)
    permeabilityRatioThreshold = 100;
end

batchDir = fileparts(mfilename('fullpath'));
rtmDir = fileparts(batchDir);
addpath(rtmDir);

runName = sprintf('exp_%03d', expIdx);
resultsDir = fullfile(batchResultsDir, runName);

cfg = struct();
cfg.runName = runName;
cfg.resultsDir = resultsDir;
cfg.outputRoot = batchResultsDir;
cfg.layoutType = params.layoutType;
cfg.characteristicLength = params.L_cm;
cfg.inletVelocity = params.u_in;
cfg.initialHydrogenConcentration = params.c_in;
cfg.permeabilityRatioThreshold = permeabilityRatioThreshold;
cfg.diffusionCoefficient = getParam(params, 'diffusionCoefficient', 1e-5);
cfg.molarVolume = getParam(params, 'molarVolume', 36.9);
cfg.rateCoefficientTST = getParam(params, 'rateCoefficientTST', 1e-4);
cfg.initialMacroscaleTimeStepSize = getParam(params, 'initialMacroscaleTimeStepSize', 0.10);
cfg.numPartitionsMicroscale = getParam(params, 'numPartitionsMicroscale', 2 * 64);
cfg.dxfResolutionX = getParam(params, 'dxfResolutionX', 200);
cfg.dxfResolutionY = getParam(params, 'dxfResolutionY', 100);
cfg.useExternalGeometry = getParam(params, 'useExternalGeometry', false);
cfg.tifPath = getParam(params, 'tifPath', "");

if isfield(params, 'Time_stepmax') && ~isempty(params.Time_stepmax)
    cfg.maximalStep = params.Time_stepmax;
end
if isfield(params, 'endTime') && ~isempty(params.endTime)
    cfg.endTime = params.endTime;
end
if isfield(params, 'targetLengthYAxis') && ~isempty(params.targetLengthYAxis)
    cfg.targetLengthYAxis = params.targetLengthYAxis;
end
if isfield(params, 'targetAspectRatio') && ~isempty(params.targetAspectRatio)
    cfg.targetAspectRatio = params.targetAspectRatio;
end
if isfield(params, 'loadExistingGeometry') && ~isempty(params.loadExistingGeometry)
    cfg.loadExistingGeometry = params.loadExistingGeometry;
end
if isfield(params, 'geometryLoadFile') && ~isempty(params.geometryLoadFile)
    cfg.geometryLoadFile = params.geometryLoadFile;
end
if isfield(params, 'geometrySaveFile') && ~isempty(params.geometrySaveFile)
    cfg.geometrySaveFile = params.geometrySaveFile;
end

% Batch defaults favor dataset generation over heavy visual output. DXF export
% stays on because the NMR/COMSOL stage consumes dxf_pore and dxf_solid.
cfg.exportEvery = getParam(params, 'exportEvery', 1);
cfg.exportDXF = getParam(params, 'exportDXF', true);
cfg.saveMainPlot = getParam(params, 'saveMainPlot', true);
cfg.saveIndividualPlots = getParam(params, 'saveIndividualPlots', false);
cfg.saveInterfaceMask = getParam(params, 'saveInterfaceMask', true);
cfg.saveRealtimePlot = getParam(params, 'saveRealtimePlot', false);
cfg.saveFigureFiles = getParam(params, 'saveFigureFiles', false);
cfg.writeExcel = getParam(params, 'writeExcel', true);
cfg.saveFinalPlot = getParam(params, 'saveFinalPlot', true);
cfg.enableNMRSimulation = getParam(params, 'enableNMRSimulation', false);

result = PNM_beauty3(cfg);
end

function value = getParam(params, fieldName, defaultValue)
if isfield(params, fieldName) && ~isempty(params.(fieldName))
    value = params.(fieldName);
else
    value = defaultValue;
end
end

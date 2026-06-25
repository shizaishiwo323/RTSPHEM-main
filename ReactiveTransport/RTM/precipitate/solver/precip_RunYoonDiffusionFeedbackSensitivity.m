function result = precip_RunYoonDiffusionFeedbackSensitivity(spec, options)
% precip_RunYoonDiffusionFeedbackSensitivity - Run Yoon n=0/2/3 short smokes.
%
% Inputs:
%   spec    - base benchmark spec.
%   options - options passed to precip_RunYoonCase1Short.
%
% Output:
%   result  - table of Case 2, Case 1, and Case 3 area outcomes.

if nargin < 1 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end
if nargin < 2
    options = struct();
end

caseNames = ["case_2"; "case_1"; "case_3"];
diffusionExponent = [0; 2; 3];
totalPrecipitatedArea_cm2 = zeros(3, 1);
totalPrecipitateMoles = zeros(3, 1);
finalMeanEffectiveDiffusivity_cm2_s = zeros(3, 1);
feedbackCoupled = true(3, 1);

for iCase = 1:numel(diffusionExponent)
    caseSpec = spec;
    caseSpec.diffusionExponent = diffusionExponent(iCase);
    caseOptions = options;
    caseOptions.coupleDiffusionFeedback = true;
    run = precip_RunYoonCase1Short(caseSpec, caseOptions);
    totalPrecipitatedArea_cm2(iCase) = ...
        run.areaTimeseries.totalPrecipitatedArea_cm2(end);
    totalPrecipitateMoles(iCase) = sum(run.state.precipitateMoles(:));
    finalMeanEffectiveDiffusivity_cm2_s(iCase) = ...
        mean(run.state.effectiveDiffusivity_cm2_s(:), 'omitnan');
end

result = table(caseNames, diffusionExponent, totalPrecipitatedArea_cm2, ...
    totalPrecipitateMoles, finalMeanEffectiveDiffusivity_cm2_s, feedbackCoupled);

manifestPath = getFieldOrDefault(options, 'outputManifestPath', '');
if ~isempty(manifestPath)
    areaTrendNonincreasing = isNonincreasing(totalPrecipitatedArea_cm2);
    massTrendDecreasing = isStrictlyDecreasing(totalPrecipitateMoles);
    manifest = struct();
    manifest.runner = 'precip_RunYoonDiffusionFeedbackSensitivity';
    manifest.caseNames = cellstr(caseNames);
    manifest.diffusionExponent = diffusionExponent(:)';
    manifest.totalPrecipitatedArea_cm2 = totalPrecipitatedArea_cm2(:)';
    manifest.totalPrecipitateMoles = totalPrecipitateMoles(:)';
    manifest.finalMeanEffectiveDiffusivity_cm2_s = ...
        finalMeanEffectiveDiffusivity_cm2_s(:)';
    manifest.feedbackCoupled = feedbackCoupled(:)';
    manifest.areaTrendNonincreasing = areaTrendNonincreasing;
    manifest.massTrendDecreasing = massTrendDecreasing;
    manifest.diffusionFeedbackAccepted = areaTrendNonincreasing && ...
        massTrendDecreasing && all(feedbackCoupled);
    writeJsonManifest(manifestPath, manifest);
end
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function tf = isNonincreasing(values)
values = double(values(:));
tolerance = 10 * eps(max(1, max(abs(values))));
tf = all(diff(values) <= tolerance);
end

function tf = isStrictlyDecreasing(values)
values = double(values(:));
tolerance = 10 * eps(max(1, max(abs(values))));
tf = all(diff(values) < -tolerance);
end

function writeJsonManifest(path, manifest)
parentDir = fileparts(path);
if ~isempty(parentDir) && ~isfolder(parentDir)
    mkdir(parentDir);
end
fid = fopen(path, 'w');
if fid < 0
    error('RTSPHEM:Precipitate:ManifestWriteFailed', ...
        'Could not open manifest for writing: %s.', path);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(manifest, 'PrettyPrint', true));
clear cleanupObj;
end

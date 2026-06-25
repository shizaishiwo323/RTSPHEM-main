function reference = precip_LoadReferenceCurves(referenceCsv, digitizationPackageDir)
% precip_LoadReferenceCurves - Load Zhang/Yoon digitized reference curves.
%
% The loader validates the repository reference-curve schema and records
% provenance metadata. It does not synthesize missing literature data.

if nargin < 1 || isempty(referenceCsv)
    referenceCsv = defaultReferenceCsv();
end
if nargin < 2
    digitizationPackageDir = [];
end
referenceCsv = char(string(referenceCsv));
if exist(referenceCsv, 'file') ~= 2
    error('RTSPHEM:Precipitate:ReferenceCurvesNotFound', ...
        'Reference curve CSV does not exist: %s', referenceCsv);
end

referenceTable = readReferenceCurves(referenceCsv);
sourceCases = unique(referenceTable(:, {'source', 'case', 'region'}), ...
    'rows', 'stable');

reference = struct();
reference.table = referenceTable;
reference.csvPath = referenceCsv;
reference.numRows = height(referenceTable);
reference.sourceCases = sourceCases;
reference.isQuantitativeBenchmark = false;
reference.note = ['Loaded approximate digitized reference curves; use ' ...
    'source screenshots, calibration files, and uncertainty records before ' ...
    'quantitative Zhang/Yoon claims.'];

if ~isempty(digitizationPackageDir)
    package = precip_LoadReferenceDigitizationPackage(digitizationPackageDir);
    if ~sameFile(package.convertedReferenceCsv, referenceCsv)
        error('RTSPHEM:Precipitate:ReferenceDigitizationPackageMismatch', ...
            ['Digitization package convertedReferenceCsv does not match ', ...
            'the loaded reference CSV.']);
    end
    reference.digitizationPackage = package;
    reference.isQuantitativeBenchmark = package.isQuantitativeBenchmark;
    if package.isQuantitativeBenchmark
        reference.note = ['Loaded reference curves with complete ', ...
            'digitization package provenance for quantitative Zhang/Yoon ', ...
            'claims.'];
    else
        reference.note = ['Loaded reference curves with incomplete ', ...
            'digitization package provenance; ', package.note];
    end
end
end

function referenceCsv = defaultReferenceCsv()
moduleRoot = fileparts(fileparts(mfilename('fullpath')));
referenceCsv = fullfile(moduleRoot, 'reference_data', ...
    'zhang_yoon_reference_curves.csv');
end

function tf = sameFile(leftPath, rightPath)
tf = strcmpi(char(string(leftPath)), char(string(rightPath)));
end

function referenceTable = readReferenceCurves(referenceCsv)
referenceTable = readtable(referenceCsv, 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');
requiredColumns = {'source', 'case', 'region', 'time_min', ...
    'precipitated_area_norm', 'precipitated_area_cm2', 'note'};
requireColumns(referenceTable, requiredColumns, ...
    'RTSPHEM:Precipitate:InvalidReferenceCurves', referenceCsv);
if height(referenceTable) == 0
    error('RTSPHEM:Precipitate:MissingReferenceCurves', ...
        'Reference curve CSV has the correct schema but no digitized rows: %s', referenceCsv);
end

referenceTable.source = string(referenceTable.source);
referenceTable.("case") = string(referenceTable.("case"));
referenceTable.region = string(referenceTable.region);
referenceTable.note = string(referenceTable.note);
referenceTable.data_role = repmat("reference", height(referenceTable), 1);

allowedSources = ["Zhang2010", "Yoon2012"];
allowedCases = ["experiment_25mM", "case_1", "case_5"];
allowedRegions = ["entire_domain", "first_pore", "first_three_pores", ...
    "zhang_upgradient", "zhang_middle", "zhang_downgradient"];
if any(~ismember(referenceTable.source, allowedSources)) || ...
        any(~ismember(referenceTable.("case"), allowedCases)) || ...
        any(~ismember(referenceTable.region, allowedRegions))
    error('RTSPHEM:Precipitate:InvalidReferenceCurves', ...
        'Reference curves contain source/case/region values outside the allowed Zhang/Yoon schema.');
end
validPairs = (referenceTable.source == "Zhang2010" & ...
    referenceTable.("case") == "experiment_25mM") | ...
    (referenceTable.source == "Yoon2012" & ...
    ismember(referenceTable.("case"), ["case_1", "case_5"]));
if any(~validPairs)
    error('RTSPHEM:Precipitate:InvalidReferenceCurves', ...
        'Reference curves contain invalid source/case pairings.');
end
validRegions = (referenceTable.source == "Zhang2010" & ...
    ismember(referenceTable.region, ...
    ["zhang_upgradient", "zhang_middle", "zhang_downgradient"])) | ...
    (referenceTable.source == "Yoon2012" & ...
    ismember(referenceTable.region, ...
    ["entire_domain", "first_pore", "first_three_pores"]));
if any(~validRegions)
    error('RTSPHEM:Precipitate:InvalidReferenceCurves', ...
        'Reference curves contain invalid source/case/region combinations.');
end
if any(~isfinite(referenceTable.time_min)) || ...
        any(~isfinite(referenceTable.precipitated_area_norm) & ...
        ~isfinite(referenceTable.precipitated_area_cm2))
    error('RTSPHEM:Precipitate:InvalidReferenceCurves', ...
        'Reference curves require finite time_min and at least one finite area value per row.');
end
referenceTable = fillMissingReferenceNormFromCm2(referenceTable);
referenceTable = referenceTable(:, {'source', 'case', 'region', 'time_min', ...
    'precipitated_area_norm', 'precipitated_area_cm2', 'note', 'data_role'});
end

function requireColumns(tableData, requiredColumns, errorId, path)
missing = setdiff(requiredColumns, tableData.Properties.VariableNames);
if ~isempty(missing)
    error(errorId, 'Missing required columns in %s: %s', path, ...
        strjoin(missing, ', '));
end
end

function referenceTable = fillMissingReferenceNormFromCm2(referenceTable)
groups = unique(referenceTable(:, {'source', 'case', 'region'}), ...
    'rows', 'stable');
for iGroup = 1:height(groups)
    mask = referenceTable.source == groups.source(iGroup) & ...
        referenceTable.("case") == groups.("case")(iGroup) & ...
        referenceTable.region == groups.region(iGroup);
    missingNorm = mask & ~isfinite(referenceTable.precipitated_area_norm) & ...
        isfinite(referenceTable.precipitated_area_cm2);
    if any(missingNorm)
        groupArea = referenceTable.precipitated_area_cm2(mask);
        scale = max(abs(groupArea), [], 'omitnan');
        if ~isfinite(scale) || scale <= 0
            error('RTSPHEM:Precipitate:InvalidReferenceCurves', ...
                'Cannot normalize reference cm2 values for %s/%s/%s.', ...
                groups.source(iGroup), groups.("case")(iGroup), ...
                groups.region(iGroup));
        end
        referenceTable.precipitated_area_norm(missingNorm) = ...
            referenceTable.precipitated_area_cm2(missingNorm) ./ scale;
        referenceTable.note(missingNorm) = referenceTable.note(missingNorm) + ...
            " | missing reference normalized values are filled from cm2 by group max absolute area";
    end
end
end

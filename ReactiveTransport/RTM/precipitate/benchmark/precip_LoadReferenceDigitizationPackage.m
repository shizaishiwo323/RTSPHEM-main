function package = precip_LoadReferenceDigitizationPackage(packageDir)
% precip_LoadReferenceDigitizationPackage - Validate reference digitization assets.
%
% A quantitative Zhang/Yoon reference package must include the source
% screenshot, WebPlotDigitizer project, axis calibration, raw export,
% conversion script, uncertainty record, and converted reference CSV. This
% loader makes that requirement machine-readable; it does not invent missing
% data or upgrade approximate rows to quantitative evidence.

if nargin < 1 || isempty(packageDir)
    error('RTSPHEM:Precipitate:MissingDigitizationPackageDir', ...
        'packageDir is required.');
end
packageDir = char(string(packageDir));
if ~isfolder(packageDir)
    error('RTSPHEM:Precipitate:DigitizationPackageNotFound', ...
        'Digitization package directory does not exist: %s', packageDir);
end

manifestPath = fullfile(packageDir, 'digitization_manifest.json');
if exist(manifestPath, 'file') ~= 2
    error('RTSPHEM:Precipitate:DigitizationManifestNotFound', ...
        'Digitization manifest does not exist: %s', manifestPath);
end

manifest = jsondecode(fileread(manifestPath));
requiredFields = {'packageName', 'sourceFigure', 'screenshotFile', ...
    'webPlotDigitizerProjectFile', 'calibrationCsv', 'rawExportCsv', ...
    'conversionScript', 'uncertaintyCsv', 'convertedReferenceCsv'};
requireFields(manifest, requiredFields, manifestPath);

assetFields = requiredFields(3:end);
assetPaths = strings(numel(assetFields), 1);
missingAssets = strings(0, 1);
for iAsset = 1:numel(assetFields)
    fieldName = assetFields{iAsset};
    assetPaths(iAsset) = string(resolvePackagePath(packageDir, ...
        manifest.(fieldName)));
    if exist(assetPaths(iAsset), 'file') ~= 2
        missingAssets(end + 1, 1) = assetPaths(iAsset); %#ok<AGROW>
    end
end

package = struct();
package.packageDir = string(packageDir);
package.manifestPath = string(manifestPath);
package.packageName = string(manifest.packageName);
package.sourceFigure = string(manifest.sourceFigure);
package.screenshotFile = assetPaths(1);
package.webPlotDigitizerProjectFile = assetPaths(2);
package.calibrationCsv = assetPaths(3);
package.rawExportCsv = assetPaths(4);
package.conversionScript = assetPaths(5);
package.uncertaintyCsv = assetPaths(6);
package.convertedReferenceCsv = assetPaths(7);
package.missingAssets = cellstr(missingAssets);

if ~isempty(missingAssets)
    package.isQuantitativeBenchmark = false;
    package.note = sprintf(['Digitization package has missing required ', ...
        'assets (%d); do not use for quantitative Zhang/Yoon claims.'], ...
        numel(missingAssets));
    return;
end

validation = validateDigitizationTables(package.calibrationCsv, ...
    package.rawExportCsv, package.uncertaintyCsv);
package.calibrationTableVerified = validation.calibrationTableVerified;
package.rawExportTableVerified = validation.rawExportTableVerified;
package.uncertaintyTableVerified = validation.uncertaintyTableVerified;
package.invalidAssets = cellstr(validation.invalidAssets);
if ~validation.accepted
    package.assetFilesVerified = false;
    package.isQuantitativeBenchmark = false;
    package.note = ['Digitization package has invalid calibration, raw ', ...
        'export, or uncertainty tables; do not use for quantitative ', ...
        'Zhang/Yoon claims.'];
    return;
end

package.reference = precip_LoadReferenceCurves(package.convertedReferenceCsv);
package.assetFilesVerified = true;
package.numReferenceRows = package.reference.numRows;
package.isQuantitativeBenchmark = true;
package.note = ['Complete digitization package with source screenshot, ', ...
    'WebPlotDigitizer project, calibration, raw export, conversion script, ', ...
    'uncertainty record, and converted reference CSV.'];
end

function requireFields(s, requiredFields, manifestPath)
missing = setdiff(requiredFields, fieldnames(s));
if ~isempty(missing)
    error('RTSPHEM:Precipitate:InvalidDigitizationManifest', ...
        'Missing required fields in %s: %s', manifestPath, ...
        strjoin(missing, ', '));
end
end

function path = resolvePackagePath(packageDir, value)
path = char(string(value));
if isempty(path)
    path = fullfile(packageDir, '__missing_empty_path__');
elseif ~isAbsolutePath(path)
    path = fullfile(packageDir, path);
end
end

function tf = isAbsolutePath(path)
tf = ~isempty(regexp(path, '^[A-Za-z]:[\\/]', 'once')) || ...
    startsWith(path, filesep);
end

function validation = validateDigitizationTables(calibrationCsv, rawExportCsv, ...
    uncertaintyCsv)
validation = struct();
validation.calibrationTableVerified = hasNumericColumns(calibrationCsv, ...
    {'pixel_x', 'pixel_y', 'data_x', 'data_y'}, 2);
validation.rawExportTableVerified = hasNumericColumns(rawExportCsv, ...
    {'x', 'y'}, 1);
validation.uncertaintyTableVerified = hasNumericColumns(uncertaintyCsv, ...
    {'uncertainty'}, 1);
validation.invalidAssets = strings(0, 1);
if ~validation.calibrationTableVerified
    validation.invalidAssets(end + 1, 1) = string(calibrationCsv);
end
if ~validation.rawExportTableVerified
    validation.invalidAssets(end + 1, 1) = string(rawExportCsv);
end
if ~validation.uncertaintyTableVerified
    validation.invalidAssets(end + 1, 1) = string(uncertaintyCsv);
end
validation.accepted = isempty(validation.invalidAssets);
end

function tf = hasNumericColumns(csvPath, requiredColumns, minRows)
tf = false;
try
    tbl = readtable(csvPath);
catch
    return;
end
if height(tbl) < minRows
    return;
end
names = string(tbl.Properties.VariableNames);
for iColumn = 1:numel(requiredColumns)
    columnName = string(requiredColumns{iColumn});
    if ~any(names == columnName)
        return;
    end
    values = tbl.(char(columnName));
    if ~isnumeric(values) || any(~isfinite(double(values(:))))
        return;
    end
    if contains(lower(columnName), "uncertainty") && any(values(:) < 0)
        return;
    end
end
tf = true;
end

function package = precip_LoadZhangYoonGeometryPackage(packageDir, spec)
% precip_LoadZhangYoonGeometryPackage - Validate Zhang/Yoon geometry assets.
%
% A quantitative geometry package must include the calibrated source image,
% axis/length calibration, processing script, substrate mask, Yoon region
% masks, and an uncertainty record. This loader validates the package and then
% delegates mask parsing to the existing geometry and region-mask loaders.

if nargin < 2
    error('RTSPHEM:Precipitate:MissingGeometrySpec', ...
        'spec is required for geometry package validation.');
end
if nargin < 1 || isempty(packageDir)
    error('RTSPHEM:Precipitate:MissingGeometryPackageDir', ...
        'packageDir is required.');
end
packageDir = char(string(packageDir));
if ~isfolder(packageDir)
    error('RTSPHEM:Precipitate:GeometryPackageNotFound', ...
        'Geometry package directory does not exist: %s', packageDir);
end

manifestPath = fullfile(packageDir, 'geometry_manifest.json');
if exist(manifestPath, 'file') ~= 2
    error('RTSPHEM:Precipitate:GeometryManifestNotFound', ...
        'Geometry manifest does not exist: %s', manifestPath);
end

manifest = jsondecode(fileread(manifestPath));
requiredFields = {'packageName', 'sourceFigure', 'sourceImageFile', ...
    'calibrationCsv', 'processingScript', 'substrateMaskFile', ...
    'regionMaskFile', 'uncertaintyCsv'};
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
package.sourceImageFile = assetPaths(1);
package.calibrationCsv = assetPaths(2);
package.processingScript = assetPaths(3);
package.substrateMaskFile = assetPaths(4);
package.regionMaskFile = assetPaths(5);
package.uncertaintyCsv = assetPaths(6);
package.missingAssets = cellstr(missingAssets);

if ~isempty(missingAssets)
    package.isQuantitativeGeometry = false;
    package.note = sprintf(['Geometry package has missing required assets ', ...
        '(%d); do not use for quantitative Zhang/Yoon geometry claims.'], ...
        numel(missingAssets));
    return;
end

validation = validateGeometryTables(package.calibrationCsv, ...
    package.uncertaintyCsv);
package.calibrationTableVerified = validation.calibrationTableVerified;
package.uncertaintyTableVerified = validation.uncertaintyTableVerified;
package.invalidAssets = cellstr(validation.invalidAssets);
if ~validation.accepted
    package.assetFilesVerified = false;
    package.isQuantitativeGeometry = false;
    package.note = ['Geometry package has invalid calibration or ', ...
        'uncertainty tables; do not use for quantitative Zhang/Yoon ', ...
        'geometry claims.'];
    return;
end

package.geometry = precip_LoadZhangYoonGeometry(package.substrateMaskFile, spec);
package.regionMasks = precip_LoadYoonRegionMasks(package.regionMaskFile, spec);
package.assetFilesVerified = true;
package.numSubstrateCells = nnz(package.geometry.substrateMask);
package.numFirstPoreCells = nnz(package.regionMasks.firstPoreMask);
package.numFirstThreePoresCells = nnz(package.regionMasks.firstThreePoresMask);
package.isQuantitativeGeometry = true;
package.note = ['Complete geometry package with calibrated source image, ', ...
    'processing script, substrate mask, region masks, and uncertainty record.'];
end

function requireFields(s, requiredFields, manifestPath)
missing = setdiff(requiredFields, fieldnames(s));
if ~isempty(missing)
    error('RTSPHEM:Precipitate:InvalidGeometryManifest', ...
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

function validation = validateGeometryTables(calibrationCsv, uncertaintyCsv)
validation = struct();
validation.calibrationTableVerified = hasNumericColumns(calibrationCsv, ...
    {'image_x_px', 'image_y_px', 'x_cm', 'y_cm'}, 2);
validation.uncertaintyTableVerified = hasNumericColumns(uncertaintyCsv, ...
    {'uncertainty_cm'}, 1);
validation.invalidAssets = strings(0, 1);
if ~validation.calibrationTableVerified
    validation.invalidAssets(end + 1, 1) = string(calibrationCsv);
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

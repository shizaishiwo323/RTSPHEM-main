function masks = precip_LoadYoonRegionMasks(maskFile, spec)
% precip_LoadYoonRegionMasks - Load Yoon area region masks from file.
%
% Supported formats:
%   .mat - variables firstPoreMask, firstThreePoresMask, optional source.
%   .csv - columns row, col, firstPoreMask, firstThreePoresMask.

if nargin < 2
    error('RTSPHEM:Precipitate:MissingRegionMaskSpec', ...
        'spec is required for region-mask size validation.');
end
if ~isfile(maskFile)
    error('RTSPHEM:Precipitate:RegionMaskFileNotFound', ...
        'Region mask file does not exist: %s.', maskFile);
end

[~, ~, ext] = fileparts(maskFile);
switch lower(ext)
    case '.mat'
        masks = loadMatMasks(maskFile);
    case '.csv'
        masks = loadCsvMasks(maskFile, spec);
    otherwise
        error('RTSPHEM:Precipitate:UnsupportedRegionMaskFile', ...
            'Unsupported region mask file extension: %s.', ext);
end

masks.firstPoreMask = validateMask(masks.firstPoreMask, spec, ...
    'firstPoreMask');
masks.firstThreePoresMask = validateMask(masks.firstThreePoresMask, spec, ...
    'firstThreePoresMask');
masks.filePath = maskFile;
if ~isfield(masks, 'source') || isempty(masks.source)
    masks.source = ['file:' maskFile];
end
end

function masks = loadMatMasks(maskFile)
data = load(maskFile);
requiredFields = {'firstPoreMask', 'firstThreePoresMask'};
for iField = 1:numel(requiredFields)
    if ~isfield(data, requiredFields{iField})
        error('RTSPHEM:Precipitate:MissingRegionMaskVariable', ...
            'MAT region mask file is missing %s.', requiredFields{iField});
    end
end
masks = struct();
masks.firstPoreMask = data.firstPoreMask;
masks.firstThreePoresMask = data.firstThreePoresMask;
if isfield(data, 'source')
    masks.source = char(string(data.source));
end
end

function masks = loadCsvMasks(maskFile, spec)
maskTable = readtable(maskFile);
requiredColumns = {'row', 'col', 'firstPoreMask', 'firstThreePoresMask'};
for iColumn = 1:numel(requiredColumns)
    if ~ismember(requiredColumns{iColumn}, maskTable.Properties.VariableNames)
        error('RTSPHEM:Precipitate:MissingRegionMaskColumn', ...
            'CSV region mask file is missing %s.', requiredColumns{iColumn});
    end
end
rows = maskTable.row;
cols = maskTable.col;
if any(rows < 1 | rows > spec.numY | cols < 1 | cols > spec.numX | ...
        rows ~= floor(rows) | cols ~= floor(cols))
    error('RTSPHEM:Precipitate:InvalidRegionMaskIndex', ...
        'CSV region mask row/col values must be valid 1-based grid indices.');
end
masks = struct();
masks.firstPoreMask = false(spec.numY, spec.numX);
masks.firstThreePoresMask = false(spec.numY, spec.numX);
linearIndex = sub2ind([spec.numY, spec.numX], rows, cols);
masks.firstPoreMask(linearIndex) = logical(maskTable.firstPoreMask);
masks.firstThreePoresMask(linearIndex) = ...
    logical(maskTable.firstThreePoresMask);
masks.source = ['file:' maskFile];
end

function mask = validateMask(mask, spec, fieldName)
if ~isequal(size(mask), [spec.numY, spec.numX])
    error('RTSPHEM:Precipitate:InvalidYoonRegionMask', ...
        '%s must have size [numY, numX].', fieldName);
end
mask = logical(mask);
end

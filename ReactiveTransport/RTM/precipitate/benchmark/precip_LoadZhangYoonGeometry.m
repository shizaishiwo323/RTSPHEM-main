function geometry = precip_LoadZhangYoonGeometry(geometryFile, spec)
% precip_LoadZhangYoonGeometry - Load static substrate geometry masks.
%
% Supported formats:
%   .mat - variable substrateMask, optional source.
%   .csv - columns row, col, substrateMask.

if nargin < 2
    error('RTSPHEM:Precipitate:MissingGeometrySpec', ...
        'spec is required for geometry size validation.');
end
if ~isfile(geometryFile)
    error('RTSPHEM:Precipitate:GeometryFileNotFound', ...
        'Geometry file does not exist: %s.', geometryFile);
end

[~, ~, ext] = fileparts(geometryFile);
switch lower(ext)
    case '.mat'
        geometry = loadMatGeometry(geometryFile);
    case '.csv'
        geometry = loadCsvGeometry(geometryFile, spec);
    otherwise
        error('RTSPHEM:Precipitate:UnsupportedGeometryFile', ...
            'Unsupported geometry file extension: %s.', ext);
end

geometry.substrateMask = validateSubstrateMask( ...
    geometry.substrateMask, spec);
geometry.filePath = geometryFile;
if ~isfield(geometry, 'source') || isempty(geometry.source)
    geometry.source = ['file:' geometryFile];
end
end

function geometry = loadMatGeometry(geometryFile)
data = load(geometryFile);
if ~isfield(data, 'substrateMask')
    error('RTSPHEM:Precipitate:MissingGeometryVariable', ...
        'MAT geometry file is missing substrateMask.');
end
geometry = struct();
geometry.substrateMask = data.substrateMask;
if isfield(data, 'source')
    geometry.source = char(string(data.source));
end
end

function geometry = loadCsvGeometry(geometryFile, spec)
geometryTable = readtable(geometryFile);
requiredColumns = {'row', 'col', 'substrateMask'};
for iColumn = 1:numel(requiredColumns)
    if ~ismember(requiredColumns{iColumn}, ...
            geometryTable.Properties.VariableNames)
        error('RTSPHEM:Precipitate:MissingGeometryColumn', ...
            'CSV geometry file is missing %s.', requiredColumns{iColumn});
    end
end
rows = geometryTable.row;
cols = geometryTable.col;
if any(rows < 1 | rows > spec.numY | cols < 1 | cols > spec.numX | ...
        rows ~= floor(rows) | cols ~= floor(cols))
    error('RTSPHEM:Precipitate:InvalidGeometryIndex', ...
        'CSV geometry row/col values must be valid 1-based grid indices.');
end
geometry = struct();
geometry.substrateMask = false(spec.numY, spec.numX);
linearIndex = sub2ind([spec.numY, spec.numX], rows, cols);
geometry.substrateMask(linearIndex) = logical(geometryTable.substrateMask);
geometry.source = ['file:' geometryFile];
end

function mask = validateSubstrateMask(mask, spec)
if ~isequal(size(mask), [spec.numY, spec.numX])
    error('RTSPHEM:Precipitate:InvalidSubstrateMask', ...
        'substrateMask must have size [numY, numX].');
end
mask = logical(mask);
end

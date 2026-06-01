function geometry = ReadLayeredDxfGeometry(dxfPath, options)
% ReadLayeredDxfGeometry - Rasterize layered DXF domain/calcite geometry.
%
% Inputs:
%   dxfPath: DXF file containing a closed domain polyline and calcite solids.
%   options.domainLayerNames: layer names for the simulation domain.
%   options.solidLayerNames: layer names for calcite solid polygons/hatches.
%   options.referenceLengthDxf/referenceLengthCm: scale calibration.
%   options.resolutionX/resolutionY: raster resolution for the signed distance.
%
% Output:
%   geometry.phiCm > 0 in calcite solid, < 0 in pore space.

if nargin < 2 || isempty(options)
    options = struct();
end
if ~isfile(dxfPath)
    error('MATLAB:ReadLayeredDxfGeometry:MissingFile', 'DXF file not found: %s', dxfPath);
end

domainLayerNames = cfggetLocal(options, 'domainLayerNames', {'domin', 'domain', 'DOMAIN'});
solidLayerNames = cfggetLocal(options, 'solidLayerNames', {'calcite'});
resolutionX = cfggetLocal(options, 'resolutionX', 400);
resolutionY = cfggetLocal(options, 'resolutionY', []);
smoothingSigmaPixels = cfggetLocal(options, 'smoothingSigmaPixels', 1.5);
outsideDomainIsSolid = logical(cfggetLocal(options, 'outsideDomainIsSolid', true));
importDirection = normalizeImportDirection(cfggetLocal(options, 'importDirection', ...
    cfggetLocal(options, 'dxfImportDirection', 'as_is')));

scaleCmPerDxfUnit = cfggetLocal(options, 'scaleCmPerDxfUnit', []);
if isempty(scaleCmPerDxfUnit)
    referenceLengthDxf = cfggetLocal(options, 'referenceLengthDxf', []);
    referenceLengthCm = cfggetLocal(options, 'referenceLengthCm', []);
    if isempty(referenceLengthDxf) || isempty(referenceLengthCm)
        error('MATLAB:ReadLayeredDxfGeometry:MissingScale', ...
            'Set either scaleCmPerDxfUnit or referenceLengthDxf/referenceLengthCm.');
    end
    scaleCmPerDxfUnit = referenceLengthCm / referenceLengthDxf;
end

entities = parseDxfClosedPolygons(dxfPath);
if isempty(entities)
    error('MATLAB:ReadLayeredDxfGeometry:NoPolygons', 'No closed polygon entities were read from: %s', dxfPath);
end

domainCandidates = filterByLayer(entities, domainLayerNames);
if isempty(domainCandidates)
    error('MATLAB:ReadLayeredDxfGeometry:NoDomain', ...
        'No domain polygon found on layers: %s', strjoin(cellstr(string(domainLayerNames)), ', '));
end
[~, domainIdx] = max(arrayfun(@polygonAbsArea, domainCandidates));
domainPolygon = domainCandidates(domainIdx);

solidCandidates = filterByLayer(entities, solidLayerNames);
if isempty(solidCandidates)
    error('MATLAB:ReadLayeredDxfGeometry:NoSolid', ...
        'No calcite solid polygon found on layers: %s', strjoin(cellstr(string(solidLayerNames)), ', '));
end

domainPolygon = transformDxfEntity(domainPolygon, importDirection);
for k = 1:numel(solidCandidates)
    solidCandidates(k) = transformDxfEntity(solidCandidates(k), importDirection);
end

domainMinX = min(domainPolygon.x);
domainMaxX = max(domainPolygon.x);
domainMinY = min(domainPolygon.y);
domainMaxY = max(domainPolygon.y);
lengthXAxis = (domainMaxX - domainMinX) * scaleCmPerDxfUnit;
lengthYAxis = (domainMaxY - domainMinY) * scaleCmPerDxfUnit;
if lengthXAxis <= 0 || lengthYAxis <= 0
    error('MATLAB:ReadLayeredDxfGeometry:InvalidDomain', 'Domain bounds are not positive in DXF: %s', dxfPath);
end

resolutionX = max(3, round(resolutionX));
if isempty(resolutionY)
    resolutionY = max(3, round(resolutionX * lengthYAxis / lengthXAxis));
else
    resolutionY = max(3, round(resolutionY));
end

dx = lengthXAxis / resolutionX;
dy = lengthYAxis / resolutionY;
xCentersCm = (0.5:1:(resolutionX - 0.5)) * dx;
yCentersCm = (0.5:1:(resolutionY - 0.5)) * dy;
[Xcm, Ycm] = meshgrid(xCentersCm, yCentersCm);

domainXcm = (domainPolygon.x - domainMinX) * scaleCmPerDxfUnit;
domainYcm = (domainPolygon.y - domainMinY) * scaleCmPerDxfUnit;
domainMask = inpolygon(Xcm, Ycm, domainXcm, domainYcm);

solidMask = false(size(domainMask));
for k = 1:numel(solidCandidates)
    solidXcm = (solidCandidates(k).x - domainMinX) * scaleCmPerDxfUnit;
    solidYcm = (solidCandidates(k).y - domainMinY) * scaleCmPerDxfUnit;
    solidMask = solidMask | inpolygon(Xcm, Ycm, solidXcm, solidYcm);
end
solidMask = solidMask & domainMask;
if outsideDomainIsSolid
    solidMask = solidMask | ~domainMask;
end

phiPixels = bwdist(~solidMask) - bwdist(solidMask);
if smoothingSigmaPixels > 0
    phiPixels = imgaussfilt(phiPixels, smoothingSigmaPixels);
end
pixelSizeCm = mean([dx, dy]);

geometry = struct();
geometry.dxfPath = string(dxfPath);
geometry.phiCm = double(phiPixels) * pixelSizeCm;
geometry.solidMask = solidMask;
geometry.domainMask = domainMask;
geometry.xCentersCm = xCentersCm;
geometry.yCentersCm = yCentersCm;
geometry.lengthXAxis = lengthXAxis;
geometry.lengthYAxis = lengthYAxis;
geometry.scaleCmPerDxfUnit = scaleCmPerDxfUnit;
geometry.domainLayer = string(domainPolygon.layer);
geometry.solidLayers = unique(string({solidCandidates.layer}));
geometry.domainBoundsRaw = [domainMinX, domainMaxX, domainMinY, domainMaxY];
geometry.resolutionX = resolutionX;
geometry.resolutionY = resolutionY;
geometry.importDirection = string(importDirection);
end

function entities = parseDxfClosedPolygons(dxfPath)
rawLines = cellstr(readlines(dxfPath));
rawLines = rawLines(~cellfun(@(s) all(ismissing(string(s))), rawLines));
if mod(numel(rawLines), 2) ~= 0
    rawLines = rawLines(1:end-1);
end
codes = str2double(strtrim(rawLines(1:2:end)));
values = strtrim(rawLines(2:2:end));

entities = struct('type', {}, 'layer', {}, 'x', {}, 'y', {});
i = 1;
while i <= numel(codes)
    if codes(i) ~= 0
        i = i + 1;
        continue;
    end

    entityType = upper(values{i});
    switch entityType
        case 'LWPOLYLINE'
            [entity, nextIdx] = parseLwpolyline(codes, values, i);
            if isClosedPolygon(entity)
                entities(end+1) = entity; %#ok<AGROW>
            end
            i = nextIdx;
        case 'POLYLINE'
            [entity, nextIdx] = parsePolyline(codes, values, i);
            if isClosedPolygon(entity)
                entities(end+1) = entity; %#ok<AGROW>
            end
            i = nextIdx;
        case 'HATCH'
            [hatchEntities, nextIdx] = parseHatch(codes, values, i);
            entities = [entities, hatchEntities]; %#ok<AGROW>
            i = nextIdx;
        otherwise
            i = i + 1;
    end
end
end

function [entity, nextIdx] = parseLwpolyline(codes, values, startIdx)
endIdx = findNextEntity(codes, startIdx + 1);
entity = emptyEntity('LWPOLYLINE');
entity.layer = readLayer(codes, values, startIdx, endIdx);
closedFlag = 0;
x = [];
y = [];
i = startIdx + 1;
while i < endIdx
    if codes(i) == 70
        closedFlag = str2double(values{i});
    elseif codes(i) == 10
        xValue = str2double(values{i});
        yValue = NaN;
        j = i + 1;
        while j < endIdx && codes(j) ~= 10 && codes(j) ~= 0
            if codes(j) == 20
                yValue = str2double(values{j});
                break;
            end
            j = j + 1;
        end
        if ~isnan(xValue) && ~isnan(yValue)
            x(end+1) = xValue; %#ok<AGROW>
            y(end+1) = yValue; %#ok<AGROW>
        end
    end
    i = i + 1;
end
if bitand(round(closedFlag), 1) == 1
    [x, y] = closePolygon(x, y);
end
entity.x = x;
entity.y = y;
nextIdx = endIdx;
end

function [entity, nextIdx] = parsePolyline(codes, values, startIdx)
entity = emptyEntity('POLYLINE');
entity.layer = readLayer(codes, values, startIdx, findNextEntity(codes, startIdx + 1));
closedFlag = 0;
x = [];
y = [];
i = startIdx + 1;
while i <= numel(codes)
    if codes(i) ~= 0
        i = i + 1;
        continue;
    end
    marker = upper(values{i});
    if strcmp(marker, 'SEQEND')
        break;
    elseif strcmp(marker, 'VERTEX')
        endIdx = findNextEntity(codes, i + 1);
        vx = NaN;
        vy = NaN;
        for j = i+1:endIdx-1
            if codes(j) == 10
                vx = str2double(values{j});
            elseif codes(j) == 20
                vy = str2double(values{j});
            elseif codes(j) == 70
                closedFlag = max(closedFlag, str2double(values{j}));
            end
        end
        if ~isnan(vx) && ~isnan(vy)
            x(end+1) = vx; %#ok<AGROW>
            y(end+1) = vy; %#ok<AGROW>
        end
        i = endIdx;
    else
        break;
    end
end
if bitand(round(closedFlag), 1) == 1 || ~isempty(x)
    [x, y] = closePolygon(x, y);
end
entity.x = x;
entity.y = y;
nextIdx = min(i + 1, numel(codes) + 1);
end

function [entities, nextIdx] = parseHatch(codes, values, startIdx)
endIdx = findNextEntity(codes, startIdx + 1);
layer = readLayer(codes, values, startIdx, endIdx);
entities = struct('type', {}, 'layer', {}, 'x', {}, 'y', {});
i = startIdx + 1;
while i < endIdx
    if codes(i) == 93
        numVertices = round(str2double(values{i}));
        x = [];
        y = [];
        i = i + 1;
        while i < endIdx && numel(x) < numVertices
            if codes(i) == 10
                xValue = str2double(values{i});
                yValue = NaN;
                j = i + 1;
                while j < endIdx && codes(j) ~= 10 && codes(j) ~= 0
                    if codes(j) == 20
                        yValue = str2double(values{j});
                        break;
                    end
                    j = j + 1;
                end
                if ~isnan(xValue) && ~isnan(yValue)
                    x(end+1) = xValue; %#ok<AGROW>
                    y(end+1) = yValue; %#ok<AGROW>
                end
            end
            i = i + 1;
        end
        [x, y] = closePolygon(x, y);
        entity = emptyEntity('HATCH');
        entity.layer = layer;
        entity.x = x;
        entity.y = y;
        if isClosedPolygon(entity)
            entities(end+1) = entity; %#ok<AGROW>
        end
    else
        i = i + 1;
    end
end
nextIdx = endIdx;
end

function nextIdx = findNextEntity(codes, startIdx)
nextIdx = numel(codes) + 1;
for k = startIdx:numel(codes)
    if codes(k) == 0
        nextIdx = k;
        return;
    end
end
end

function layer = readLayer(codes, values, startIdx, endIdx)
layer = '0';
for k = startIdx+1:endIdx-1
    if codes(k) == 8
        layer = values{k};
        return;
    end
end
end

function out = filterByLayer(entities, layerNames)
wanted = lower(string(layerNames));
keep = false(size(entities));
for k = 1:numel(entities)
    keep(k) = any(lower(string(entities(k).layer)) == wanted);
end
out = entities(keep);
end

function areaValue = polygonAbsArea(entity)
if numel(entity.x) < 3
    areaValue = 0;
else
    areaValue = abs(polyarea(entity.x, entity.y));
end
end

function tf = isClosedPolygon(entity)
tf = numel(entity.x) >= 4 && entity.x(1) == entity.x(end) && entity.y(1) == entity.y(end);
end

function entity = emptyEntity(typeName)
entity = struct('type', typeName, 'layer', '0', 'x', [], 'y', []);
end

function [x, y] = closePolygon(x, y)
if isempty(x)
    return;
end
if x(1) ~= x(end) || y(1) ~= y(end)
    x(end+1) = x(1);
    y(end+1) = y(1);
end
end

function direction = normalizeImportDirection(direction)
direction = lower(strtrim(char(direction)));
direction = strrep(direction, '-', '_');
direction = strrep(direction, ' ', '_');
switch direction
    case {'', 'asis', 'as_is', 'none', 'original'}
        direction = 'as_is';
    case {'rotate90_cw', 'rotate90cw', 'rotate_90_cw', 'rot90cw'}
        direction = 'rotate90_cw';
    case {'rotate90_ccw', 'rotate90ccw', 'rotate_90_ccw', 'rot90ccw'}
        direction = 'rotate90_ccw';
    case {'rotate180', 'rotate_180', 'rot180'}
        direction = 'rotate180';
    case {'flipx', 'flip_x', 'mirror_x'}
        direction = 'flip_x';
    case {'flipy', 'flip_y', 'mirror_y'}
        direction = 'flip_y';
    otherwise
        error('MATLAB:ReadLayeredDxfGeometry:UnknownImportDirection', ...
            'Unknown DXF importDirection "%s".', direction);
end
end

function entity = transformDxfEntity(entity, direction)
x = entity.x;
y = entity.y;
switch direction
    case 'as_is'
        xNew = x;
        yNew = y;
    case 'rotate90_cw'
        xNew = y;
        yNew = -x;
    case 'rotate90_ccw'
        xNew = -y;
        yNew = x;
    case 'rotate180'
        xNew = -x;
        yNew = -y;
    case 'flip_x'
        xNew = -x;
        yNew = y;
    case 'flip_y'
        xNew = x;
        yNew = -y;
    otherwise
        error('MATLAB:ReadLayeredDxfGeometry:UnknownImportDirection', ...
            'Unknown DXF importDirection "%s".', direction);
end
entity.x = xNew;
entity.y = yNew;
end

function value = cfggetLocal(config, fieldName, defaultValue)
if isstruct(config) && isfield(config, fieldName) && ~isempty(config.(fieldName))
    value = config.(fieldName);
else
    value = defaultValue;
end
end

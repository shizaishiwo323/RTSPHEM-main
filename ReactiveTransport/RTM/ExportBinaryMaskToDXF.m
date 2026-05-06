function ExportBinaryMaskToDXF(X, Y, mask, filename, layerName)
    if nargin < 5 || isempty(layerName)
        layerName = '0';
    end
    if all(mask(:) == 0) || all(mask(:) == 1)
        return;
    end

    mask = mask ~= 0;
    maskPad = padarray(mask, [1, 1], 0, 'both');
    boundaries = bwboundaries(maskPad, 8, 'holes');
    if isempty(boundaries)
        return;
    end

    xCenters = X(1, :);
    yCenters = Y(:, 1);
    dx = mean(diff(xCenters));
    dy = mean(diff(yCenters));

    xEdges = (xCenters(1) - dx/2) + dx * (0:(numel(xCenters)+1));
    yEdges = (yCenters(1) - dy/2) + dy * (0:(numel(yCenters)+1));

    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open %s for writing.', filename);
    end

    fprintf(fid, '0\nSECTION\n2\nHEADER\n0\nENDSEC\n');
    fprintf(fid, '0\nSECTION\n2\nENTITIES\n');

    for k = 1:numel(boundaries)
        outline = boundaries{k};
        rowIdx = outline(:, 1) - 1;
        colIdx = outline(:, 2) - 1;

        if any(colIdx < 0 | colIdx > numel(xEdges)-1 | ...
               rowIdx < 0 | rowIdx > numel(yEdges)-1)
            continue;
        end

        x = xEdges(colIdx + 1).';
        y = yEdges(rowIdx + 1).';

        if numel(x) < 3
            continue;
        end
        if x(1) ~= x(end) || y(1) ~= y(end)
            x(end+1) = x(1);
            y(end+1) = y(1);
        end

        fprintf(fid, '0\nPOLYLINE\n8\n%s\n66\n1\n70\n1\n', layerName);
        for i = 1:numel(x)
            fprintf(fid, '0\nVERTEX\n8\n%s\n10\n%.6f\n20\n%.6f\n30\n0.0\n', layerName, x(i), y(i));
        end
        fprintf(fid, '0\nSEQEND\n');
    end

    fprintf(fid, '0\nENDSEC\n0\nEOF\n');
    fclose(fid);
end
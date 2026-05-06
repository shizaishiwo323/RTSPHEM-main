function ExportDomainRectangleToDXF(lengthX, lengthY, filename, layerName)
    if nargin < 4 || isempty(layerName)
        layerName = '0';
    end
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open %s for writing.', filename);
    end
    fprintf(fid, '0\nSECTION\n2\nENTITIES\n');
    fprintf(fid, '0\nPOLYLINE\n8\n%s\n66\n1\n70\n1\n', layerName);
    verts = [0, 0; lengthX, 0; lengthX, lengthY; 0, lengthY; 0, 0];
    for i = 1:size(verts, 1)
        fprintf(fid, '0\nVERTEX\n8\n%s\n10\n%.6f\n20\n%.6f\n30\n0.0\n', layerName, verts(i, 1), verts(i, 2));
    end
    fprintf(fid, '0\nSEQEND\n0\nENDSEC\n0\nEOF\n');
    fclose(fid);
end
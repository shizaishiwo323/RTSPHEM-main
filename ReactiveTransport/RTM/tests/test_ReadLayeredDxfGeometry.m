function tests = test_ReadLayeredDxfGeometry
tests = functiontests(localfunctions);
end

function testReadsLayeredDxfAndScalesToCentimeters(testCase)
dxfPath = [tempname, '.dxf'];
cleanup = onCleanup(@() deleteIfExists(dxfPath));
writeFixtureDxf(dxfPath);

options = struct();
options.domainLayerNames = {'domin'};
options.solidLayerNames = {'calcite'};
options.referenceLengthDxf = 1500;
options.referenceLengthCm = 0.15;
options.resolutionX = 80;
options.resolutionY = 50;
options.smoothingSigmaPixels = 0;

geometry = ReadLayeredDxfGeometry(dxfPath, options);

verifyEqual(testCase, geometry.lengthXAxis, 0.15, 'AbsTol', 1e-12);
verifyEqual(testCase, geometry.lengthYAxis, 0.10, 'AbsTol', 1e-12);
verifySize(testCase, geometry.phiCm, [50, 80]);
verifyClass(testCase, geometry.phiCm, 'double');
verifyTrue(testCase, geometry.solidMask(25, 35));
verifyFalse(testCase, geometry.solidMask(5, 5));
verifyGreaterThan(testCase, geometry.phiCm(25, 35), 0);
verifyLessThan(testCase, geometry.phiCm(5, 5), 0);
end

function testCanRotateDxfImportDirection(testCase)
dxfPath = [tempname, '.dxf'];
cleanup = onCleanup(@() deleteIfExists(dxfPath));
writeFixtureDxf(dxfPath);

options = struct();
options.domainLayerNames = {'domin'};
options.solidLayerNames = {'calcite'};
options.referenceLengthDxf = 1500;
options.referenceLengthCm = 0.15;
options.resolutionX = 50;
options.resolutionY = 80;
options.smoothingSigmaPixels = 0;
options.importDirection = 'rotate90_ccw';

geometry = ReadLayeredDxfGeometry(dxfPath, options);

verifyEqual(testCase, geometry.lengthXAxis, 0.10, 'AbsTol', 1e-12);
verifyEqual(testCase, geometry.lengthYAxis, 0.15, 'AbsTol', 1e-12);
verifyTrue(testCase, any(geometry.solidMask(:)));
verifyEqual(testCase, geometry.importDirection, "rotate90_ccw");
end

function writeFixtureDxf(path)
fid = fopen(path, 'w');
assert(fid ~= -1, 'Could not create DXF fixture.');
fprintf(fid, '0\nSECTION\n2\nENTITIES\n');
writeLwpolyline(fid, 'domin', [0 1500 1500 0], [0 0 1000 1000]);
writeHatchPolygon(fid, 'calcite', [500 900 900 500], [300 300 700 700]);
fprintf(fid, '0\nENDSEC\n0\nEOF\n');
fclose(fid);
end

function writeLwpolyline(fid, layerName, x, y)
fprintf(fid, '0\nLWPOLYLINE\n8\n%s\n90\n%d\n70\n1\n', layerName, numel(x));
for i = 1:numel(x)
    fprintf(fid, '10\n%.12g\n20\n%.12g\n', x(i), y(i));
end
end

function writeHatchPolygon(fid, layerName, x, y)
fprintf(fid, '0\nHATCH\n8\n%s\n91\n1\n92\n7\n72\n0\n73\n1\n93\n%d\n', layerName, numel(x));
for i = 1:numel(x)
    fprintf(fid, '10\n%.12g\n20\n%.12g\n', x(i), y(i));
end
end

function deleteIfExists(path)
if exist(path, 'file')
    delete(path);
end
end

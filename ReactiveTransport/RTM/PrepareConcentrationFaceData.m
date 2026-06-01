function faceData = PrepareConcentrationFaceData(concentrationData, triangles, levelSetData)
% PrepareConcentrationFaceData - Use P0 triangle concentration and mask solids.

faceData = concentrationData(:);
if isempty(levelSetData)
    return;
end

levelSetOnTriangles = levelSetData(triangles);
poreTriangles = all(levelSetOnTriangles < 0, 2);
faceData(~poreTriangles) = NaN;
end

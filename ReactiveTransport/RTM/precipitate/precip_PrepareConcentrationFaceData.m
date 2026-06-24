function faceData = precip_PrepareConcentrationFaceData(concentrationData, triangles, levelSetData)
% precip_PrepareConcentrationFaceData - Use P0 triangle values and mask solids.

faceData = concentrationData(:);
if isempty(levelSetData)
    return;
end

levelSetOnTriangles = levelSetData(triangles);
poreTriangles = all(levelSetOnTriangles < 0, 2);
faceData(~poreTriangles) = NaN;
end

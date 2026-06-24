function ratePerArea = ComputeLegacyEquivalentTstTriangleRate( ...
    triangleHydrogenMolCm3, interfaceAreaCm2, gridV0T, vertexTriMatrix, ...
    numVertices, rateCoefficientTST)
% Compute triangle rates using the legacy RTSPHEM TST vertex discretization.
%
% The legacy moving-boundary model first assigns each vertex the minimum H+
% value over adjacent triangles, then uses r = k * c_H * 1000 at the
% interface. This helper maps that same vertex rate back to triangles for
% PHREEQC prescribed-reaction moles.

triangleHydrogenMolCm3 = max(triangleHydrogenMolCm3(:), 0);
interfaceAreaCm2 = max(interfaceAreaCm2(:), 0);
vertexRate = zeros(numVertices, 1);
for iVertex = 1:numVertices
    triIndices = vertexTriMatrix{1, iVertex};
    values = triangleHydrogenMolCm3(triIndices);
    values = values(isfinite(values));
    if isempty(values)
        vertexRate(iVertex) = 0;
    else
        vertexRate(iVertex) = min(values) * 1000 * rateCoefficientTST;
    end
end
ratePerArea = mean(vertexRate(gridV0T), 2);
ratePerArea(interfaceAreaCm2 <= 0) = 0;
ratePerArea(~isfinite(ratePerArea)) = 0;
end

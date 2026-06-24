function step = ComputeTstMatchStepDissolution( ...
    triangleRateDriverMolCm3, interfaceAreaCm2, calciteMoles, timeStepSize, ...
    gridV0T, vertexTriMatrix, numVertices, rateCoefficientTST, molarVolumeCm3Mol)
% Compute one explicit TST-match dissolution step from a single rate driver.
%
% The same vertex-minimum rate driver is used for moving-boundary normal
% speed and for the PHREEQC prescribed CaCO3 moles.  For legacy TST the
% driver is c_H; for PHREEQC-TST it is a_H = gamma_H*c_H from pH.

triangleRateDriverMolCm3 = max(triangleRateDriverMolCm3(:), 0);
interfaceAreaCm2 = max(interfaceAreaCm2(:), 0);
calciteMoles = max(calciteMoles(:), 0);
timeStepSize = max(timeStepSize, 0);

vertexRateMolCm2S = zeros(numVertices, 1);
for iVertex = 1:numVertices
    triIndices = vertexTriMatrix{1, iVertex};
    values = triangleRateDriverMolCm3(triIndices);
    values = values(isfinite(values));
    if isempty(values)
        vertexRateMolCm2S(iVertex) = 0;
    else
        vertexRateMolCm2S(iVertex) = min(values) * 1000 * rateCoefficientTST;
    end
end

ratePerArea = mean(vertexRateMolCm2S(gridV0T), 2);
ratePerArea(interfaceAreaCm2 <= 0) = 0;
ratePerArea(~isfinite(ratePerArea)) = 0;

prescribedMoles = ratePerArea(:) .* interfaceAreaCm2(:) .* timeStepSize;
prescribedMoles = min(max(prescribedMoles, 0), calciteMoles(:));
prescribedMoles(~isfinite(prescribedMoles)) = 0;

normalSpeed = molarVolumeCm3Mol .* vertexRateMolCm2S(:);
normalSpeed(~isfinite(normalSpeed)) = 0;

step = struct();
step.ratePerArea_mol_cm2_s = ratePerArea(:);
step.prescribedMoles = prescribedMoles(:);
step.normalSpeed_cm_s = normalSpeed(:);
step.vertexRate_mol_cm2_s = vertexRateMolCm2S(:);
end

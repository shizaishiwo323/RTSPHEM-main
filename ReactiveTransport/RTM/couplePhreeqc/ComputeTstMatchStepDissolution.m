function step = ComputeTstMatchStepDissolution( ...
    triangleRateDriverMolCm3, interfaceAreaCm2, calciteMoles, timeStepSize, ...
    gridV0T, vertexTriMatrix, numVertices, rateCoefficientTST, molarVolumeCm3Mol)
% Compute one explicit TST-match dissolution step from a single rate driver.
%
% The same vertex-minimum rate driver is used to propose CaCO3 dissolution.
% The realized, inventory-limited moles are then the single source of truth
% for both the PHREEQC prescribed reaction and moving-boundary normal speed.

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

triangleVertexRates = reshape(vertexRateMolCm2S(gridV0T(:)), size(gridV0T));
ratePerArea = mean(triangleVertexRates, 2);
ratePerArea(interfaceAreaCm2 <= 0) = 0;
ratePerArea(~isfinite(ratePerArea)) = 0;

candidateMoles = ratePerArea(:) .* interfaceAreaCm2(:) .* timeStepSize;
candidateMoles = max(candidateMoles, 0);
candidateMoles(~isfinite(candidateMoles)) = 0;

realizedMoles = min(candidateMoles(:), calciteMoles(:));
realizedMoles(~isfinite(realizedMoles)) = 0;
inventoryLimited = candidateMoles(:) > realizedMoles(:);

realizedRatePerArea = zeros(size(ratePerArea(:)));
active = interfaceAreaCm2(:) > 0 & timeStepSize > 0;
realizedRatePerArea(active) = realizedMoles(active) ./ ...
    interfaceAreaCm2(active) ./ timeStepSize;
realizedRatePerArea(~isfinite(realizedRatePerArea)) = 0;

normalSpeed = zeros(numVertices, 1);
for iVertex = 1:numVertices
    triIndices = vertexTriMatrix{1, iVertex};
    values = realizedRatePerArea(triIndices);
    values = values(isfinite(values));
    if isempty(values)
        normalSpeed(iVertex) = 0;
    else
        normalSpeed(iVertex) = mean(values) * molarVolumeCm3Mol;
    end
end
normalSpeed(~isfinite(normalSpeed)) = 0;

step = struct();
step.candidateRatePerArea_mol_cm2_s = ratePerArea(:);
step.ratePerArea_mol_cm2_s = realizedRatePerArea(:);
step.candidateMoles = candidateMoles(:);
step.realizedMoles = realizedMoles(:);
step.prescribedMoles = realizedMoles(:);
step.inventoryLimited = inventoryLimited(:);
step.reactantLimited = false(size(realizedMoles(:)));
step.normalSpeed_cm_s = normalSpeed(:);
step.vertexRate_mol_cm2_s = vertexRateMolCm2S(:);
step.realizedRatePerArea_mol_cm2_s = realizedRatePerArea(:);
end

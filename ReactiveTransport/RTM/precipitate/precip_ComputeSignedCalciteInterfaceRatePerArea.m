function ratePerArea = precip_ComputeSignedCalciteInterfaceRatePerArea(result, interfaceAreaCm2, timeStepSize)
% precip_ComputeSignedCalciteInterfaceRatePerArea - Dissolution-positive rate.

if isfield(result, 'calciteDeltaMoles') && ~isempty(result.calciteDeltaMoles)
    deltaMoles = result.calciteDeltaMoles(:);
elseif isfield(result, 'calciteSignedRate_mol_s') && ~isempty(result.calciteSignedRate_mol_s)
    deltaMoles = result.calciteSignedRate_mol_s(:) .* max(timeStepSize, eps);
else
    deltaMoles = zeros(size(interfaceAreaCm2(:)));
end

interfaceAreaCm2 = interfaceAreaCm2(:);
ratePerArea = zeros(size(deltaMoles));
active = interfaceAreaCm2 > 0 & timeStepSize > 0;
ratePerArea(active) = -deltaMoles(active) ./ max(timeStepSize, eps) ./ interfaceAreaCm2(active);
ratePerArea(~isfinite(ratePerArea)) = 0;
end

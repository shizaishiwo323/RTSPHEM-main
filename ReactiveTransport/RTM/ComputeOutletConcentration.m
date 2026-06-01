function outletConcentration = ComputeOutletConcentration(concentrationData, outletTriangleMask)
% ComputeOutletConcentration - Mean concentration in outlet-adjacent cells.

mask = outletTriangleMask(:);
if ~any(mask)
    outletConcentration = NaN;
    return;
end

outletConcentration = mean(concentrationData(mask), 'omitnan');
end

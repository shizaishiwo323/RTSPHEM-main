function samples = precip_BuildYoonMixingSeries(spec, fractionInletA)
% precip_BuildYoonMixingSeries - Build zero-dimensional split-inlet mixtures.
%
% Inputs:
%   spec            - benchmark spec.
%   fractionInletA  - fraction of CaCl2 inlet in each mixture.
%
% Output:
%   samples         - component totals and endpoint pH metadata.

fractionInletA = fractionInletA(:);
samples = struct();
samples.fractionInletA = fractionInletA;
samples.fixedPH = nan(size(fractionInletA));
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    samples.(fieldName) = fractionInletA .* spec.inletA.(fieldName) + ...
        (1 - fractionInletA) .* spec.inletB.(fieldName);
end
samples.fixedPH(abs(fractionInletA - 1) < 10 * eps) = spec.inletA.pH;
samples.fixedPH(abs(fractionInletA) < 10 * eps) = spec.inletB.pH;
end

function rate = precip_YoonVateriteRate(chem, spec)
% precip_YoonVateriteRate - Evaluate the signed Yoon/Chou Vaterite rate.
%
% Inputs:
%   chem - struct with aH, aH2CO3, and omegaVaterite arrays.
%   spec - benchmark spec with yoonRate constants.
%
% Output:
%   rate - mol/cm2/s; positive values precipitate, negative values dissolve.

if nargin < 2 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end
requiredFields = {'aH', 'aH2CO3', 'omegaVaterite'};
for iField = 1:numel(requiredFields)
    if ~isfield(chem, requiredFields{iField})
        error('RTSPHEM:Precipitate:MissingChemistryField', ...
            'Chemistry field %s is required.', requiredFields{iField});
    end
end

rateConstants = spec.yoonRate;
prefactor = rateConstants.k1 .* chem.aH + ...
    rateConstants.k2 .* chem.aH2CO3 + rateConstants.k3;
rate = -prefactor .* (1 - chem.omegaVaterite);
if isfield(spec, 'dissolutionFactor') && isfinite(spec.dissolutionFactor)
    dissolutionMask = rate < 0;
    rate(dissolutionMask) = spec.dissolutionFactor .* rate(dissolutionMask);
end
end

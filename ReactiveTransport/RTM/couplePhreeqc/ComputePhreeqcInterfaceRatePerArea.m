function ratePerArea = ComputePhreeqcInterfaceRatePerArea(result, interfaceAreaCm2, timeStepSize, options, legacyHMolCm3)
% ComputePhreeqcInterfaceRatePerArea - Rate used to move calcite interface.
%
% For external_tst_phreeqc, the transported H+ field defines the first-order
% interface law r_area = c_H * 1000 * kTST. The caller prescribes the
% corresponding CaCO3 amount to PHREEQC and reuses this rate field for the
% moving boundary.

if nargin < 5
    legacyHMolCm3 = [];
end

rateLaw = lower(strrep(strtrim(char(getOption(options, 'rateLaw', ...
    getOption(options, 'phreeqcRateLaw', 'database_calcite')))), '-', '_'));

switch rateLaw
    case {'tst_match', 'calcite_tst_match', 'phreeqc_tst_match', ...
            'external_tst_phreeqc', 'legacy_phreeqc_tst_match'}
        rateCoefficientTST = getOption(options, 'rateCoefficientTST', ...
            getOption(options, 'phreeqcTstRateCoefficient', 1e-4));
        if isempty(legacyHMolCm3)
            if isfield(result, 'h_activity_mol_cm3') && ~isempty(result.h_activity_mol_cm3)
                legacyHMolCm3 = result.h_activity_mol_cm3;
            else
                legacyHMolCm3 = result.h_mol_cm3;
            end
        end
        ratePerArea = max(legacyHMolCm3(:), 0) * 1000 * rateCoefficientTST;
        active = interfaceAreaCm2(:) > 0;
        ratePerArea(~active) = 0;
    otherwise
        if isfield(result, 'calciteDissolvedMoles')
            dissolvedMoles = result.calciteDissolvedMoles;
        else
            dissolvedMoles = max(result.calciteRate_mol_s * timeStepSize, 0);
        end
        ratePerArea = zeros(size(dissolvedMoles));
        active = interfaceAreaCm2(:) > 0 & timeStepSize > 0;
        ratePerArea(active) = dissolvedMoles(active) ./ max(timeStepSize, eps) ./ interfaceAreaCm2(active);
end

ratePerArea(~isfinite(ratePerArea)) = 0;
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

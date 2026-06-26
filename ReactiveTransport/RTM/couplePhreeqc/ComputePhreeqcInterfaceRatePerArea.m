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
        rawMotionRate = rawCalciteKineticRate(result, timeStepSize);
        if ~isempty(rawMotionRate)
            currentMolarVolume = getOption(options, 'molarVolume', ...
                getOption(options, 'molarVolume_cm3_mol', 36.9));
            legacySpaceScale = getOption(options, 'phreeqcLegacySpeedScaleFactor', 20);
            legacyCalciteMolarVolume = getOption(options, ...
                'phreeqcLegacyCalciteMolarVolumeForSpeed', 39.63e-3);
            ratePerArea = rawMotionRate .* legacySpaceScale .* ...
                legacyCalciteMolarVolume ./ max(currentMolarVolume, eps);
            ratePerArea(interfaceAreaCm2(:) <= 0) = 0;
            ratePerArea(~isfinite(ratePerArea)) = 0;
            return;
        end
        interfaceRateMode = normalizeRateMode(getOption(options, ...
            'phreeqcInterfaceRateMode', getOption(options, ...
            'interfaceRateMode', 'per_area_from_dissolved_moles')));
        if strcmp(interfaceRateMode, 'kin_delta_over_kin_time')
            if isfield(result, 'calciteKinDeltaRate_mol_s') && ...
                    ~isempty(result.calciteKinDeltaRate_mol_s)
                ratePerArea = max(result.calciteKinDeltaRate_mol_s(:), 0);
            elseif isfield(result, 'calciteRate_mol_s') && ~isempty(result.calciteRate_mol_s)
                ratePerArea = max(result.calciteRate_mol_s(:), 0);
            elseif isfield(result, 'calciteDissolvedMoles')
                ratePerArea = max(result.calciteDissolvedMoles(:), 0) ./ ...
                    max(timeStepSize, eps);
            else
                ratePerArea = zeros(size(interfaceAreaCm2(:)));
            end
            ratePerArea(interfaceAreaCm2(:) <= 0) = 0;
            ratePerArea(~isfinite(ratePerArea)) = 0;
            return;
        end
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

function rawRate = rawCalciteKineticRate(result, timeStepSize)
rawRate = [];
if isfield(result, 'calciteRawKineticDeltaMoles') && ...
        ~isempty(result.calciteRawKineticDeltaMoles)
    rawRate = max(-result.calciteRawKineticDeltaMoles(:), 0) ./ ...
        max(timeStepSize, eps);
elseif isfield(result, 'calciteRawKinDeltaRate_mol_s') && ...
        ~isempty(result.calciteRawKinDeltaRate_mol_s)
    rawRate = max(result.calciteRawKinDeltaRate_mol_s(:), 0);
end
if ~isempty(rawRate) && ~any(rawRate > 0 & isfinite(rawRate))
    rawRate = [];
end
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

function mode = normalizeRateMode(mode)
mode = lower(strrep(strtrim(char(mode)), '-', '_'));
switch mode
    case {'kin_delta_over_kin_time', 'kin_delta_rate', ...
            'phreeqc_kin_delta_rate', 'phreeqc_cell_rate', ...
            'direct_kin_delta'}
        mode = 'kin_delta_over_kin_time';
    otherwise
        mode = 'per_area_from_dissolved_moles';
end
end

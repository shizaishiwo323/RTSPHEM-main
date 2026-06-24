function [upperBounds, speciesNames] = precip_ComputePhreeqcTransportUpperBounds(config, maxFactor)
% precip_ComputePhreeqcTransportUpperBounds
% Compute per-species transport caps before passing concentrations to PHREEQC.

speciesNames = {'H', 'Ca_total', 'C_total', 'Na_total', 'Cl_total'};
if nargin < 2 || isempty(maxFactor)
    maxFactor = Inf;
end
if ~isscalar(maxFactor) || ~isnumeric(maxFactor) || isnan(maxFactor) || maxFactor <= 0
    error('RTSPHEM:Precipitate:InvalidTransportMaxFactor', ...
        'phreeqcTransportMaxFactor must be a positive numeric scalar or Inf.');
end
if isinf(maxFactor)
    upperBounds = Inf(numel(speciesNames), 1);
    return;
end

baseMaxima = [
    maxConcentration(config, 'initialHydrogenConcentration', '', 'H_total')
    maxConcentration(config, 'initialCalciumConcentration', 'inletCalciumConcentration', 'Ca_total')
    maxConcentration(config, 'initialCarbonConcentration', 'inletCarbonConcentration', 'C_total')
    maxConcentration(config, 'initialSodiumConcentration', 'inletSodiumConcentration', 'Na_total')
    maxConcentration(config, 'initialChlorideConcentration', 'inletChlorideConcentration', 'Cl_total')
    ];
upperBounds = max(baseMaxima, 0) .* maxFactor;
end

function value = maxConcentration(config, initialField, scalarInletField, splitField)
values = zeros(1, 4);
values(1) = scalarField(config, initialField);
if ~isempty(scalarInletField)
    values(2) = scalarField(config, scalarInletField);
end
values(3) = nestedScalarField(config, 'inletA', splitField);
values(4) = nestedScalarField(config, 'inletB', splitField);
values = values(isfinite(values));
if isempty(values)
    value = 0;
else
    value = max(values);
end
end

function value = scalarField(config, fieldName)
if isstruct(config) && isfield(config, fieldName) && ~isempty(config.(fieldName)) && ...
        isnumeric(config.(fieldName)) && isscalar(config.(fieldName))
    value = double(config.(fieldName));
else
    value = 0;
end
end

function value = nestedScalarField(config, parentField, fieldName)
if isstruct(config) && isfield(config, parentField) && isstruct(config.(parentField)) && ...
        isfield(config.(parentField), fieldName) && ~isempty(config.(parentField).(fieldName)) && ...
        isnumeric(config.(parentField).(fieldName)) && isscalar(config.(parentField).(fieldName))
    value = double(config.(parentField).(fieldName));
else
    value = 0;
end
end

function ValidatePhreeqcTransportBasis(state, options)
%VALIDATEPHREEQCTRANSPORTBASIS Reject derived species as transported components.
%
% PHREEQC-coupled transport must carry conserved totals. Free H+, pH,
% carbonate aqueous species, activities, and saturation indices are derived
% chemistry outputs and must not be written into state.component_names.

if nargin < 2 || isempty(options)
    options = struct();
end
rtm.state.ValidateState(state);

componentNames = string(state.component_names(:));
canonicalNames = lower(strrep(strrep(strtrim(componentNames), " ", ""), "_", ""));
forbidden = forbiddenTransportComponents(options);
[isForbidden, forbiddenIndex] = ismember(canonicalNames, forbidden.canonical);
if any(isForbidden)
    first = find(isForbidden, 1, 'first');
    label = forbidden.display(forbiddenIndex(first));
    error('RTSPHEM:Chemistry:ForbiddenTransportComponent', ...
        ['PHREEQC transport state contains derived component "%s". ', ...
        'Transport conserved totals instead and keep pH/free H/carbonate species ', ...
        'as PHREEQC-derived outputs.'], char(label));
end
end

function forbidden = forbiddenTransportComponents(options)
baseDisplay = ["H+", "free_H", "pH", "HCO3-", "HCO3", ...
    "CO3-2", "CO3--", "CO3", "OH-", "m_H+", "a_H+", ...
    "SI_Calcite", "calciteSI"];
if isstruct(options) && isfield(options, 'derived') && ~isempty(options.derived)
    baseDisplay = unique([baseDisplay(:); string(options.derived(:))], 'stable');
end
forbidden = struct();
forbidden.display = baseDisplay(:);
forbidden.canonical = lower(strrep(strrep(strtrim(forbidden.display), " ", ""), "_", ""));
end

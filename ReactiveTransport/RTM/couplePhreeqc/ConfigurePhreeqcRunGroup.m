function cfg = ConfigurePhreeqcRunGroup(cfg, groupName, runStamp)
% ConfigurePhreeqcRunGroup - Apply named PHREEQC chemistry group settings.
%
% Supported groups:
%   phreeqc_database_calcite : PHREEQC database Calcite RATES.
%   external_tst_phreeqc     : transport-defined first-order TST calcite
%                              amount prescribed to PHREEQC for speciation.
%
% Deprecated aliases:
%   phreeqc_tst_match        : alias for external_tst_phreeqc.

if nargin < 1 || isempty(cfg)
    cfg = struct();
end
if nargin < 2 || isempty(groupName)
    groupName = 'phreeqc_database_calcite';
end
if nargin < 3
    runStamp = '';
end

[canonicalGroup, legacyAlias] = normalizeGroupName(groupName);
if ~isempty(legacyAlias)
    warning('RTSPHEM:Phreeqc:DeprecatedRunGroup', ...
        ['PHREEQC run group "%s" is deprecated. Use "%s" for ' ...
        'explicit external TST rate + PHREEQC equilibrium closure.'], ...
        char(groupName), canonicalGroup);
end
cfg.reactionModel = 'phreeqc';
cfg.phreeqcRunGroup = canonicalGroup;
if isfield(cfg, 'legacyPhreeqcRunGroupAlias')
    cfg = rmfield(cfg, 'legacyPhreeqcRunGroupAlias');
end

switch canonicalGroup
    case 'phreeqc_database_calcite'
        cfg.phreeqcRateLaw = 'database_calcite';
        cfg.chemistryMode = 'phreeqc_kinetics';
        cfg.chemistrySemantics = 'PHREEQC phreeqc_rates.dat Calcite kinetics';
    case 'external_tst_phreeqc'
        cfg.phreeqcRateLaw = 'tst_match';
        cfg.chemistryMode = 'external_tst_phreeqc';
        cfg.chemistrySemantics = 'explicit external TST rate + PHREEQC equilibrium closure';
        if ~isempty(legacyAlias)
            cfg.legacyPhreeqcRunGroupAlias = legacyAlias;
        end
        if ~isfield(cfg, 'rateCoefficientTST') || isempty(cfg.rateCoefficientTST)
            cfg.rateCoefficientTST = getFieldOrDefault(cfg, 'phreeqcTstRateCoefficient', 1e-4);
        end
        cfg.phreeqcTstRateCoefficient = cfg.rateCoefficientTST;
        cfg.phreeqcMaxSpecificSurfaceArea = Inf;
    otherwise
        error('RTSPHEM:Phreeqc:UnknownRunGroup', ...
            'Unknown PHREEQC run group: %s.', char(groupName));
end

if ~isempty(runStamp)
    cfg.runName = sprintf('%s_%s', canonicalGroup, char(runStamp));
end
end

function [canonicalGroup, legacyAlias] = normalizeGroupName(groupName)
legacyAlias = '';
normalized = lower(strrep(strtrim(char(groupName)), '-', '_'));
switch normalized
    case {'phreeqc_database_calcite', 'database_calcite', 'database', 'calcite'}
        canonicalGroup = 'phreeqc_database_calcite';
    case {'external_tst_phreeqc', 'external_tst', 'tst_phreeqc'}
        canonicalGroup = 'external_tst_phreeqc';
    case {'phreeqc_tst_match', 'tst_match', 'calcite_tst_match', ...
            'legacy_phreeqc_tst_match'}
        canonicalGroup = 'external_tst_phreeqc';
        legacyAlias = normalized;
    otherwise
        error('RTSPHEM:Phreeqc:UnknownRunGroup', ...
            'Unknown PHREEQC run group: %s.', char(groupName));
end
end

function value = getFieldOrDefault(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function cfg = ConfigurePhreeqcRunGroup(cfg, groupName, runStamp)
% ConfigurePhreeqcRunGroup - Apply named PHREEQC chemistry group settings.
%
% Supported groups:
%   phreeqc_database_calcite : PHREEQC database Calcite RATES.
%   phreeqc_tst_match        : transport-defined first-order TST calcite
%                              amount prescribed to PHREEQC for speciation.

if nargin < 1 || isempty(cfg)
    cfg = struct();
end
if nargin < 2 || isempty(groupName)
    groupName = 'phreeqc_database_calcite';
end
if nargin < 3
    runStamp = '';
end

canonicalGroup = normalizeGroupName(groupName);
cfg.reactionModel = 'phreeqc';
cfg.phreeqcRunGroup = canonicalGroup;

switch canonicalGroup
    case 'phreeqc_database_calcite'
        cfg.phreeqcRateLaw = 'database_calcite';
    case 'phreeqc_tst_match'
        cfg.phreeqcRateLaw = 'tst_match';
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

function canonicalGroup = normalizeGroupName(groupName)
normalized = lower(strrep(strtrim(char(groupName)), '-', '_'));
switch normalized
    case {'phreeqc_database_calcite', 'database_calcite', 'database', 'calcite'}
        canonicalGroup = 'phreeqc_database_calcite';
    case {'phreeqc_tst_match', 'tst_match', 'calcite_tst_match'}
        canonicalGroup = 'phreeqc_tst_match';
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

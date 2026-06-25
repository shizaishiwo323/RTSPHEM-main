function databasePath = ResolvePhreeqcDatabasePath(rtmDir, preferredDatabaseName, options)
% ResolvePhreeqcDatabasePath - Locate the PHREEQC database for RTM runs.
%
% The benchmark and single-run entry points should be reproducible from the
% RTM tree.  Therefore the bundled runtime database under
% ReactiveTransport/RTM/phreeqc/database is preferred over user Downloads.

if nargin < 1 || isempty(rtmDir)
    rtmDir = fileparts(fileparts(mfilename('fullpath')));
end
if nargin < 2 || isempty(preferredDatabaseName)
    preferredDatabaseName = 'phreeqc-m.dat';
end
if nargin < 3 || isempty(options)
    options = struct();
end
databasePolicy = lower(strrep(strtrim(char(getOption(options, ...
    'databasePolicy', 'allow_fallback'))), '-', '_'));

runtimeDatabaseDir = fullfile(rtmDir, 'phreeqc', 'database');
exactLocalPath = fullfile(runtimeDatabaseDir, char(preferredDatabaseName));
if strcmpi(databasePolicy, 'exact_local')
    if exist(exactLocalPath, 'file') == 2
        databasePath = string(exactLocalPath);
        return;
    end
    error('RTSPHEM:Phreeqc:MissingExactLocalDatabase', ...
        ['Cannot locate the exact local PHREEQC database required by this ', ...
        'run: %s. exact_local mode does not allow fallback databases.'], ...
        exactLocalPath);
elseif ~strcmpi(databasePolicy, 'allow_fallback')
    error('RTSPHEM:Phreeqc:UnknownDatabasePolicy', ...
        'Unknown PHREEQC database policy: %s.', databasePolicy);
end

candidates = {
    exactLocalPath
    fullfile(runtimeDatabaseDir, 'phreeqc-m.dat')
    fullfile(runtimeDatabaseDir, 'phreeqc.dat')
    fullfile(rtmDir, 'couplePhreeqc', 'phreeqc-m.dat')
    char(preferredDatabaseName)
    'C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat'
    };

databasePath = "";
seen = strings(0, 1);
for iCandidate = 1:numel(candidates)
    candidate = string(candidates{iCandidate});
    if any(strcmpi(seen, candidate))
        continue;
    end
    seen(end + 1, 1) = candidate; %#ok<AGROW>
    if exist(char(candidate), 'file') == 2
        databasePath = candidate;
        return;
    end
end

error('RTSPHEM:Phreeqc:MissingDatabase', ...
    ['Cannot locate a PHREEQC database. Expected the local runtime copy at ', ...
    '%s, or an installed IPhreeqcCOM database.'], ...
        fullfile(runtimeDatabaseDir, char(preferredDatabaseName)));
end

function value = getOption(options, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
else
    value = defaultValue;
end
end

function databasePath = ResolvePhreeqcDatabasePath(rtmDir, preferredDatabaseName)
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

runtimeDatabaseDir = fullfile(rtmDir, 'phreeqc', 'database');
candidates = {
    char(preferredDatabaseName)
    fullfile(runtimeDatabaseDir, char(preferredDatabaseName))
    fullfile(runtimeDatabaseDir, 'phreeqc-m.dat')
    fullfile(runtimeDatabaseDir, 'phreeqc.dat')
    fullfile(rtmDir, 'couplePhreeqc', 'phreeqc-m.dat')
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

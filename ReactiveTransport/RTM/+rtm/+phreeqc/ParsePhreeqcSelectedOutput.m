function result = ParsePhreeqcSelectedOutput(rawOutput, expectedRows)
%PARSEPHREEQCSELECTEDOUTPUT Package entry for PHREEQC selected output parsing.

if nargin < 2
    expectedRows = [];
end
ensureLegacyCouplePhreeqcPath();
if isempty(expectedRows)
    result = ParsePhreeqcSelectedOutput(rawOutput);
else
    result = ParsePhreeqcSelectedOutput(rawOutput, expectedRows);
end
end

function ensureLegacyCouplePhreeqcPath()
rtmDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
helperDir = fullfile(rtmDir, 'couplePhreeqc');
if exist(helperDir, 'dir') == 7
    addpath(helperDir);
end
end

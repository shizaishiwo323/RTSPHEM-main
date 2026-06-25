function text = BuildCalcitePhreeqcInput(state, options)
%BUILDCALCITEPHREEQCINPUT Package entry for calcite PHREEQC input text.

if nargin < 2 || isempty(options)
    options = struct();
end
ensureLegacyCouplePhreeqcPath();
text = BuildCalcitePhreeqcInput(state, options);
end

function ensureLegacyCouplePhreeqcPath()
rtmDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
helperDir = fullfile(rtmDir, 'couplePhreeqc');
if exist(helperDir, 'dir') == 7
    addpath(helperDir);
end
end

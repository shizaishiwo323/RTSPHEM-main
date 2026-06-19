function variableObj = ExtendHyPHMVariableToStepper(variableObj, fillMode)
% ExtendHyPHMVariableToStepper expands a HyPHM Variable after its Stepper
% time grid has been extended.

if nargin < 2 || isempty(fillMode)
    fillMode = 'empty';
end
fillMode = char(fillMode);

requiredContentLength = variableObj.stepper.numsteps + 1;
currentContentLength = numel(variableObj.content);
if currentContentLength < requiredContentLength
    for iEntry = (currentContentLength + 1):requiredContentLength
        switch fillMode
            case 'copy_previous'
                if currentContentLength >= 1
                    variableObj.content{iEntry, 1} = variableObj.content{currentContentLength, 1};
                else
                    variableObj.content{iEntry, 1} = [];
                end
            case 'empty'
                variableObj.content{iEntry, 1} = [];
            otherwise
                error('RTM:InvalidVariableExtensionMode', ...
                    'Unsupported fillMode: %s', fillMode);
        end
    end
end

if isfieldOrProp(variableObj, 'grids')
    requiredGridLength = variableObj.stepper.numsteps;
    currentGridLength = numel(variableObj.grids);
    if currentGridLength < requiredGridLength
        for iEntry = (currentGridLength + 1):requiredGridLength
            if currentGridLength >= 1
                variableObj.grids{iEntry, 1} = variableObj.grids{currentGridLength, 1};
            else
                variableObj.grids{iEntry, 1} = [];
            end
        end
    end
end
end

function tf = isfieldOrProp(value, name)
if isstruct(value)
    tf = isfield(value, name);
else
    tf = isprop(value, name);
end
end

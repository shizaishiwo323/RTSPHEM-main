function summary = precip_ExtendTransportProblemsToStepper(transportProblems, extensionFcn)
% precip_ExtendTransportProblemsToStepper - Extend transport variables after appending steps.

if nargin < 2 || isempty(extensionFcn)
    extensionFcn = @ExtendHyPHMVariableToStepper;
end
if ~iscell(transportProblems)
    transportProblems = {transportProblems};
end

variableFields = {'Q', 'U', 'A', 'B', 'C', 'D', 'E', 'F', ...
    'uD', 'gF', 'uR', 'rob', 'L'};

numVariablesExtended = 0;
for iProblem = 1:numel(transportProblems)
    transportProblem = transportProblems{iProblem};
    for iField = 1:numel(variableFields)
        fieldName = variableFields{iField};
        [hasValue, value] = getTransportVariable(transportProblem, fieldName);
        if hasValue && ~isempty(value)
            extensionFcn(value, 'copy_previous');
            numVariablesExtended = numVariablesExtended + 1;
        end
    end
end

summary = struct( ...
    'numProblems', numel(transportProblems), ...
    'numVariablesExtended', numVariablesExtended);
end

function [hasValue, value] = getTransportVariable(transportProblem, fieldName)
hasValue = false;
value = [];
if isstruct(transportProblem)
    if isfield(transportProblem, fieldName)
        hasValue = true;
        value = transportProblem.(fieldName);
    end
elseif isobject(transportProblem)
    if isprop(transportProblem, fieldName)
        hasValue = true;
        value = transportProblem.(fieldName);
    end
end
end

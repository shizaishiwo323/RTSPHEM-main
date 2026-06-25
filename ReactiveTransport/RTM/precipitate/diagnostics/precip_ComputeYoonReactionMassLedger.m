function ledger = precip_ComputeYoonReactionMassLedger(initialState, currentState, spec)
% precip_ComputeYoonReactionMassLedger - Audit Yoon aqueous/solid inventories.
%
% Conserved quantities for the Vaterite reaction are:
%   Ca_total + precipitate
%   C_total + precipitate
%   Alkalinity + 2 * precipitate
% Na and Cl are reaction-inert and are reported as aqueous inventories.

initialSolid = totalPrecipitateMoles(initialState);
currentSolid = totalPrecipitateMoles(currentState);

before = struct();
before.totalCa_mol = aqueousMass(initialState, 'Ca_total', spec) + initialSolid;
before.totalC_mol = aqueousMass(initialState, 'C_total', spec) + initialSolid;
before.totalAlkalinityEq = aqueousMass(initialState, 'Alkalinity', spec) + ...
    2 * initialSolid;
before.Na_total_mol = aqueousMass(initialState, 'Na_total', spec);
before.Cl_total_mol = aqueousMass(initialState, 'Cl_total', spec);
before.precipitateMoles = initialSolid;

after = struct();
after.totalCa_mol = aqueousMass(currentState, 'Ca_total', spec) + currentSolid;
after.totalC_mol = aqueousMass(currentState, 'C_total', spec) + currentSolid;
after.totalAlkalinityEq = aqueousMass(currentState, 'Alkalinity', spec) + ...
    2 * currentSolid;
after.Na_total_mol = aqueousMass(currentState, 'Na_total', spec);
after.Cl_total_mol = aqueousMass(currentState, 'Cl_total', spec);
after.precipitateMoles = currentSolid;

ledger = struct();
ledger.before = before;
ledger.after = after;
ledger.delta = subtractStruct(after, before);
ledger.relative = relativeStruct(ledger.delta, before);
conservedRelative = [ledger.relative.totalCa_mol, ...
    ledger.relative.totalC_mol, ledger.relative.totalAlkalinityEq, ...
    ledger.relative.Na_total_mol, ledger.relative.Cl_total_mol];
ledger.accepted = max(abs(conservedRelative)) <= 1e-10;
end

function value = aqueousMass(state, fieldName, spec)
if isfield(state, 'componentMoles') && ...
        isfield(state.componentMoles, fieldName) && ...
        ~isempty(state.componentMoles.(fieldName))
    value = sum(state.componentMoles.(fieldName)(:));
    return;
end
state = precip_RefreshYoonComponentMolesFromAqueous(state, spec);
value = sum(state.componentMoles.(fieldName)(:));
end

function value = totalPrecipitateMoles(state)
if isfield(state, 'precipitateMoles') && ~isempty(state.precipitateMoles)
    value = sum(state.precipitateMoles(:));
else
    value = 0;
end
end

function out = subtractStruct(a, b)
fields = fieldnames(a);
out = struct();
for iField = 1:numel(fields)
    fieldName = fields{iField};
    out.(fieldName) = a.(fieldName) - b.(fieldName);
end
end

function out = relativeStruct(delta, before)
fields = fieldnames(delta);
out = struct();
for iField = 1:numel(fields)
    fieldName = fields{iField};
    out.(fieldName) = delta.(fieldName) / max(abs(before.(fieldName)), eps);
end
end

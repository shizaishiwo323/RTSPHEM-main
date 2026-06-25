function [updated, ledger] = precip_AdvanceMicrocontinuumTransport2D(state, spec, dt, options)
% precip_AdvanceMicrocontinuumTransport2D - Advance state component moles.
%
% The conservative variable is state.componentMoles. Aqueous concentrations
% are recovered from the current fluid volume after the inventory update.

if nargin < 4
    options = struct();
end
if dt < 0
    error('RTSPHEM:Precipitate:InvalidTransportStep', ...
        'Transport time step must be nonnegative.');
end

state = precip_RefreshYoonComponentMolesFromAqueous(state, spec);
updated = state;
ledger = initializeLedger(spec);

for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    [fluxes, scalarLedger] = precip_AssembleComponentFaceFluxes( ...
        state, fieldName, spec, options);
    moles = state.componentMoles.(fieldName);
    moles = applyFaceFluxDivergence(moles, fluxes.fluxX_mol_s, ...
        fluxes.fluxY_mol_s, dt);
    updated.componentMoles.(fieldName) = moles;
    scalarLedger.massFinal = sum(moles(:));
    scalarLedger.massChange = scalarLedger.massFinal - ...
        scalarLedger.massInitial;
    scalarLedger.boundaryClosureError = scalarLedger.massChange - ...
        scalarLedger.boundaryNetFlux * dt;
    ledger.massInitial.(fieldName) = scalarLedger.massInitial;
    ledger.massFinal.(fieldName) = scalarLedger.massFinal;
    ledger.massChange.(fieldName) = scalarLedger.massChange;
    ledger.boundaryNetFlux.(fieldName) = scalarLedger.boundaryNetFlux;
    ledger.boundaryClosureError.(fieldName) = ...
        scalarLedger.boundaryClosureError;
end
updated = precip_RefreshYoonAqueousFromComponentMoles(updated, spec);
ledger.maxBoundaryClosureError = max(abs(struct2array( ...
    ledger.boundaryClosureError)));
end

function moles = applyFaceFluxDivergence(moles, fluxX, fluxY, dt)
moles = moles + dt .* (fluxX(:, 1:end-1) - fluxX(:, 2:end));
moles = moles + dt .* (fluxY(1:end-1, :) - fluxY(2:end, :));
end

function ledger = initializeLedger(spec)
ledger = struct();
ledger.massInitial = struct();
ledger.massFinal = struct();
ledger.massChange = struct();
ledger.boundaryNetFlux = struct();
ledger.boundaryClosureError = struct();
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    ledger.massInitial.(fieldName) = NaN;
    ledger.massFinal.(fieldName) = NaN;
    ledger.massChange.(fieldName) = NaN;
    ledger.boundaryNetFlux.(fieldName) = NaN;
    ledger.boundaryClosureError.(fieldName) = NaN;
end
ledger.maxBoundaryClosureError = NaN;
end

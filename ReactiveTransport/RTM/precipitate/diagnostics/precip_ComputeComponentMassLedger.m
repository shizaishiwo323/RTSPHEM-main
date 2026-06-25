function ledger = precip_ComputeComponentMassLedger(beforeComponents, afterComponents, spec)
% precip_ComputeComponentMassLedger - Compare conservative component inventory.
%
% Inputs:
%   beforeComponents - struct of component concentration fields, mol/cm3.
%   afterComponents  - struct of component concentration fields, mol/cm3.
%   spec             - benchmark spec with component names and cell volume.
%
% Output:
%   ledger           - absolute and relative inventory differences by component.

ledger = struct();
ledger.before = struct();
ledger.after = struct();
ledger.delta = struct();
ledger.relative = struct();

for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    beforeMass = sum(beforeComponents.(fieldName)(:)) * spec.cellVolume_cm3;
    afterMass = sum(afterComponents.(fieldName)(:)) * spec.cellVolume_cm3;
    deltaMass = afterMass - beforeMass;
    ledger.before.(fieldName) = beforeMass;
    ledger.after.(fieldName) = afterMass;
    ledger.delta.(fieldName) = deltaMass;
    ledger.relative.(fieldName) = deltaMass / max(abs(beforeMass), eps);
end
end

function updated = precip_UpdateVateriteVolumeFraction(state, rate, areaCm2, dt, spec)
% precip_UpdateVateriteVolumeFraction - Apply signed Vaterite precipitation.
%
% Inputs:
%   state   - Yoon state with components and Vm.
%   rate    - signed mol/cm2/s rate; positive precipitates.
%   areaCm2 - reactive area per cell, cm2.
%   dt      - reaction time step, s.
%   spec    - benchmark spec.
%
% Output:
%   updated - state with Vm, components, fluid fraction, and precipitateMoles.
%             blockedMask is refreshed later by precip_UpdateFlowMask so that
%             topology changes remain observable by the flow step.

if dt < 0
    error('RTSPHEM:Precipitate:InvalidReactionStep', ...
        'Reaction time step must be nonnegative.');
end

updated = state;
componentVolumeCm3 = spec.cellVolume_cm3 .* ones(size(state.Vm));
deltaMoles = rate .* areaCm2 .* dt;
deltaMoles = limitDeltaMoles(deltaMoles, state, componentVolumeCm3, spec);
deltaVm = spec.vateriteMolarVolume_cm3_mol .* deltaMoles ./ spec.cellVolume_cm3;

updated.Vm = min(max(state.Vm + deltaVm, 0), 1);
acceptedDeltaMoles = (updated.Vm - state.Vm) .* spec.cellVolume_cm3 ./ ...
    spec.vateriteMolarVolume_cm3_mol;
deltaConcentration = acceptedDeltaMoles ./ componentVolumeCm3;

updated.components.Ca_total = state.components.Ca_total - deltaConcentration;
updated.components.C_total = state.components.C_total - deltaConcentration;
updated.components.Alkalinity = state.components.Alkalinity - 2 .* deltaConcentration;
updated.fluidVolumeFraction = double(~state.substrateMask) .* (1 - updated.Vm);
updated.precipitateMoles = updated.Vm .* spec.cellVolume_cm3 ./ ...
    spec.vateriteMolarVolume_cm3_mol;
end

function limitedDeltaMoles = limitDeltaMoles(deltaMoles, state, componentVolumeCm3, spec)
limitedDeltaMoles = deltaMoles;
precipitationMask = deltaMoles > 0;
dissolutionMask = deltaMoles < 0;

if any(precipitationMask(:))
    caAvailable = max(state.components.Ca_total, 0) .* componentVolumeCm3;
    cAvailable = max(state.components.C_total, 0) .* componentVolumeCm3;
    alkAvailable = max(state.components.Alkalinity, 0) .* componentVolumeCm3 ./ 2;
    maxByAqueous = min(cat(3, caAvailable, cAvailable, alkAvailable), [], 3);
    maxByVm = maxVmStepMoles(spec, state.Vm);
    maxPrecipitation = min(maxByAqueous, maxByVm);
    limitedDeltaMoles(precipitationMask) = min(deltaMoles(precipitationMask), ...
        maxPrecipitation(precipitationMask));
end

if any(dissolutionMask(:))
    precipitateMoles = state.Vm .* spec.cellVolume_cm3 ./ ...
        spec.vateriteMolarVolume_cm3_mol;
    maxByVm = maxVmStepMoles(spec, state.Vm);
    maxDissolution = min(precipitateMoles, maxByVm);
    limitedDeltaMoles(dissolutionMask) = max(deltaMoles(dissolutionMask), ...
        -maxDissolution(dissolutionMask));
end
end

function maxMoles = maxVmStepMoles(spec, vm)
if isfield(spec, 'maxVmChangePerStep') && isfinite(spec.maxVmChangePerStep)
    maxVmChange = max(spec.maxVmChangePerStep, 0);
else
    maxVmChange = Inf;
end
maxMoles = maxVmChange .* spec.cellVolume_cm3 ./ spec.vateriteMolarVolume_cm3_mol .* ...
    ones(size(vm));
end

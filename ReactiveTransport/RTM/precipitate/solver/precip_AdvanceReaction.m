function updated = precip_AdvanceReaction(state, spec, dt)
% precip_AdvanceReaction - Fixed-geometry Yoon reaction step.
%
% Inputs:
%   state - Yoon micro-continuum state with component fields and Vm.
%   spec  - benchmark spec.
%   dt    - reaction time step, s.
%
% Output:
%   updated - state after speciation, Vaterite reaction, Vm update, and Deff.

if dt < 0
    error('RTSPHEM:Precipitate:InvalidReactionStep', ...
        'Reaction time step must be nonnegative.');
end

samples = stateComponentsToSamples(state, spec);
[chem, chemistryMetadata] = precip_SpeciateYoonComponents(samples, spec);
rate = reshape(precip_YoonVateriteRate(chem, spec), size(state.Vm));
area = precip_ComputeYoonReactiveArea(state, spec);

updated = precip_UpdateVateriteVolumeFraction(state, rate, area, dt, spec);
updated.effectiveDiffusivity_cm2_s = precip_UpdateEffectiveDiffusivity(updated, spec);
updated.chemistry = reshapeChemistry(chem, size(state.Vm));
updated.chemistryBackend = chemistryMetadata.backend;
updated.chemistryBackendMetadata = chemistryMetadata;
updated.diagnostics = struct();
updated.diagnostics.reactiveArea_cm2 = area;
updated.diagnostics.rate_mol_cm2_s = rate;
updated.diagnostics.dt_s = dt;
end

function samples = stateComponentsToSamples(state, spec)
samples = struct();
samples.fixedPH = nan(numel(state.Vm), 1);
for iComponent = 1:numel(spec.componentNames)
    fieldName = spec.componentNames{iComponent};
    samples.(fieldName) = state.components.(fieldName)(:);
end
end

function chemOut = reshapeChemistry(chem, targetSize)
chemOut = struct();
fields = fieldnames(chem);
for iField = 1:numel(fields)
    fieldName = fields{iField};
    chemOut.(fieldName) = reshape(chem.(fieldName), targetSize);
end
end

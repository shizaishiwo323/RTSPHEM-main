function result = MockPhreeqcBatch(batchState, ~)
%MOCKPHREEQCBATCH Deterministic PHREEQC batch stand-in for benchmark tests.

dissolved = batchState.prescribed_calcite_dissolved_moles(:);
if isfield(batchState, 'reaction_water_volume_cm3') && ...
        ~isempty(batchState.reaction_water_volume_cm3)
    waterVolume = batchState.reaction_water_volume_cm3(:);
else
    waterVolume = batchState.water_volume_cm3(:);
end
deltaConcentration = zeros(size(dissolved));
active = waterVolume > 0;
deltaConcentration(active) = dissolved(active) ./ waterVolume(active);

result = struct();
result.ca_total_mol_cm3 = batchState.ca_total_mol_cm3(:) + deltaConcentration;
result.c_total_mol_cm3 = batchState.c_total_mol_cm3(:) + deltaConcentration;
result.na_total_mol_cm3 = batchState.na_total_mol_cm3(:);
result.cl_total_mol_cm3 = batchState.cl_total_mol_cm3(:);
result.calciteDissolvedMoles = dissolved;
result.pH = repmat(7, numel(dissolved), 1);
result.calciteSI = zeros(numel(dissolved), 1);
end

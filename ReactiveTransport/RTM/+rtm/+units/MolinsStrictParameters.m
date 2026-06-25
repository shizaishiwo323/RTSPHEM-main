function params = MolinsStrictParameters()
%MOLINSSTRICTPARAMETERS Canonical cgs parameters for strict Molins tests.

params = struct();
params.domain_length_cm = 0.1;
params.domain_height_cm = 0.05;
params.circle_radius_cm = 0.01;
params.circle_center_cm = [0.05, 0.025];
params.inlet_velocity_cm_s = 0.12;
params.diffusion_coefficient_cm2_s = 1e-5;
params.h_inlet_mol_cm3 = 1.255e-6;
params.rate_constant_cm_s = 0.1;
params.calcite_molar_volume_cm3_mol = rtm.units.CalciteMolarVolumeCm3Mol();
params.units = rtm.units.NormalizeRtmUnits();
end

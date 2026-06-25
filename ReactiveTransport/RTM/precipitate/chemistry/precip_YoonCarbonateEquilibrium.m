function chem = precip_YoonCarbonateEquilibrium(samples, spec)
% precip_YoonCarbonateEquilibrium - Minimal carbonate equilibrium for Yoon smokes.
%
% Inputs:
%   samples - component totals from precip_BuildYoonMixingSeries, mol/cm3.
%   spec    - benchmark spec with Vaterite constants.
%
% Output:
%   chem    - pH, activities, Omega and SI for Vaterite.

if nargin < 2 || isempty(spec)
    spec = precip_ZhangYoonBenchmarkSpec();
end

numSamples = numel(samples.Ca_total);
pH = zeros(numSamples, 1);
aCa = zeros(numSamples, 1);
aCO3 = zeros(numSamples, 1);
aH = zeros(numSamples, 1);
aH2CO3 = zeros(numSamples, 1);

for iSample = 1:numSamples
    if isfield(samples, 'fixedPH') && isfinite(samples.fixedPH(iSample))
        pH(iSample) = samples.fixedPH(iSample);
    else
        pH(iSample) = solvePHFromAlkalinity(samples.C_total(iSample), ...
            samples.Alkalinity(iSample));
    end
    species = carbonateSpecies(samples.C_total(iSample), pH(iSample));
    aCa(iSample) = max(samples.Ca_total(iSample) * 1000, 0);
    aCO3(iSample) = max(species.co3_mol_l, 0);
    aH(iSample) = 10^(-pH(iSample));
    aH2CO3(iSample) = max(species.h2co3_mol_l, 0);
end

omega = (aCa .* aCO3) ./ spec.vateriteKsp;
chem = struct();
chem.pH = pH;
chem.aH = aH;
chem.aCa = aCa;
chem.aCO3 = aCO3;
chem.aH2CO3 = aH2CO3;
chem.omegaVaterite = omega;
chem.saturationIndexVaterite = log10(max(omega, realmin));
end

function pH = solvePHFromAlkalinity(carbonMolCm3, alkalinityMolCm3)
if carbonMolCm3 <= 0 && abs(alkalinityMolCm3) < 1e-18
    pH = 7;
    return;
end

targetAlk = alkalinityMolCm3 * 1000;
lo = 0;
hi = 14;
for iter = 1:100
    mid = 0.5 * (lo + hi);
    residual = alkalinityAtPH(carbonMolCm3, mid) - targetAlk;
    if residual > 0
        hi = mid;
    else
        lo = mid;
    end
end
pH = 0.5 * (lo + hi);
end

function alkMolL = alkalinityAtPH(carbonMolCm3, pH)
species = carbonateSpecies(carbonMolCm3, pH);
h = 10^(-pH);
kw = 1e-14;
alkMolL = species.hco3_mol_l + 2 * species.co3_mol_l + kw / h - h;
end

function species = carbonateSpecies(carbonMolCm3, pH)
carbonMolL = max(carbonMolCm3 * 1000, 0);
h = 10^(-pH);
k1 = 10^(-6.35);
k2 = 10^(-10.33);
denom = h.^2 + k1 .* h + k1 .* k2;
species = struct();
species.h2co3_mol_l = carbonMolL .* h.^2 ./ denom;
species.hco3_mol_l = carbonMolL .* k1 .* h ./ denom;
species.co3_mol_l = carbonMolL .* k1 .* k2 ./ denom;
end

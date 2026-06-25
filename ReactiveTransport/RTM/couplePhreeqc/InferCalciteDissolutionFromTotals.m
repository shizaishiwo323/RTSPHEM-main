function result = InferCalciteDissolutionFromTotals(result, state, timeStepSize)
%INFERCALCITEDISSOLUTIONFROMTOTALS Compatibility wrapper for diagnostics.
%
% Production code must not infer mineral deltas from aqueous totals. This
% legacy entry now delegates to DiagnosticInferCalciteDissolutionFromTotals
% and preserves PHREEQC-reported calciteDissolvedMoles/calciteDeltaMoles.

if nargin < 3
    timeStepSize = [];
end
result = DiagnosticInferCalciteDissolutionFromTotals(result, state, timeStepSize);
end

# Zhang/Yoon Vaterite Benchmark Sources

This directory records source-bounded constants for the Zhang 2010 / Yoon
2012 calcium-carbonate precipitation benchmark. The MATLAB specification in
`precip_ZhangYoonBenchmarkSpec.m` is the executable source of truth for the
new Yoon micro-continuum path.

## Current Scope

- Mineral phase: Vaterite.
- Solid state: per-cell precipitate volume fraction `Vm`.
- Conservative transport components: `Ca_total`, `C_total`, `Na_total`,
  `Cl_total`, and `Alkalinity`.
- Split inlets: 25 mM CaCl2 at pH 6.1 and 25 mM Na2CO3 at pH 10.9.
- Base transport values: Darcy velocity 1.25 cm/min and aqueous diffusion
  coefficient 9.0e-10 m2/s.
- Yoon feedback defaults: `D_eff = D * (1 - Vm)^2` and blocked cells at
  `Vm >= 0.6`.
- Default Vaterite rate constants are locked to the Yoon/Chou values converted
  to centimeter units: `k1 = 8.9e-5`, `k2 = 5.01e-8`, and
  `k3 = 6.6e-11 mol_cm-2_s-1`. Vaterite solubility uses
  `Ksp = 1.832e-8 mol2_L-2` and `1.832e-2 mol2_m-6`; molar volume is
  `37.47 cm3/mol`.
  Fixed-geometry manifests record `yoonRateK1`, `yoonRateK2`, `yoonRateK3`,
  `yoonRateUnits`, `yoonRateSource`, `yoonRateSourceDoi`,
  `yoonRateSourceEquation`, `yoonRateMineralPhase`, and
  `yoonRateSourceValuesVerified`; readiness audits reject locked-rate claims
  that omit these finite constants, use the wrong units, omit source
  provenance, target a non-Vaterite phase, or have not verified the source
  values.

## Boundary Between Model Paths

The existing signed Calcite/level-set code remains available for diagnostic
surface-growth smokes. It should not be used to claim a strict Zhang/Yoon
Vaterite reproduction.

The signed Calcite PHREEQC builder no longer defaults to a pH 7 floor for
positive strong-base H+ inputs; its default `minHForPHMolL` is 1e-12 so the
25 mM Na2CO3 inlet pH 10.9 can be represented when that legacy path is run.
This is a bug fix for the diagnostic path, not a change in benchmark scope.
The same legacy driver now initializes the H+ transport variable from the
configured `initialHydrogenConcentration` instead of hard-coding zero. Free H+
transport remains outside the recommended Yoon micro-continuum path.

The new Yoon path starts from static substrate geometry plus dynamic `Vm`.
PHREEQC may later be used for aqueous speciation, but the benchmark reaction
law and solid update are owned by the Yoon/Vaterite micro-continuum solver.
Alkalinity is a signed conservative component in this path: PHREEQC speciation
inputs preserve negative values, and the Vaterite precipitation limiter uses
Ca/C inventories rather than nonnegative Alkalinity availability.

The current solver foundation includes top/bottom reactive area, static
substrate-wall area, exposed `Vm = 1` precipitate-face area without internal
filled-cell double counting, signed Vaterite rate evaluation, stoichiometric
Ca/C/alkalinity updates, dynamic `Vm`, `D_eff = D * (1 - Vm)^n` calculation,
and a fixed-geometry reaction step that chains these operators. A chemistry
backend dispatcher now makes the aqueous
speciation role explicit: `yoon_equilibrium` is the current local carbonate
equilibrium smoke backend, while `iphreeqc_speciation` is an injection interface
for a caller-supplied `spec.iphreeqcSpeciationFcn`. The supplied
`precip_IPhreeqcSpeciation` implementation builds zero-dimensional PHREEQC
solutions from the Yoon conservative components and parses pH, activities, and
Vaterite SI through a persistent `rtm.phreeqc.PhreeqcSession`. It defines
`PHASES Vaterite` from `spec.vateriteKsp` inside the PHREEQC input because the
bundled `phreeqc.dat` lacks that phase. This path is unit-tested with a mock
PHREEQC engine, and a local three-point IPhreeqcCOM smoke has run with a
supersaturated mixed sample. A 101-point local IPhreeqcCOM mixing-series smoke
also runs through `precip_RunYoonSpeciationMixingSeries`; PHREEQC currently
peaks near `fractionInletA = 0.45` with `SI_Vaterite` about 2.51, whereas the
local Yoon equilibrium peak is at 0.5. Quantitative PHREEQC/Yoon equilibrium
cross-validation remains pending because this discrepancy has not yet been
accepted or calibrated. The mixing-series runner can now write a JSON chemistry
manifest with pH/SI differences, thresholds, PHREEQC backend provenance,
`numSamples`, `mixingFractionStep`, `mixingFractionsCover101`, and
`chemistrySpeciationAccepted`; readiness audits require that flag plus
`hasIphreeqc`, `isQuantitativeAcceptance`, full `0:0.01:1` coverage, and finite
pH/SI differences within their finite thresholds before quantitative claims.
The chemistry gate also
requires PHREEQC inlet pH endpoints to preserve the expected pH 6.1/10.9 values
within the recorded inlet-pH tolerance and requires both Yoon and PHREEQC
Vaterite peaks to occur at an internal mixing fraction, not at an inlet endpoint. A short
fixed-geometry Case 1 smoke driver is present for center-band and
mass-diagnostic testing. A reaction mass ledger records closure of
`Ca_total + precipitate`, `C_total + precipitate`, and
`Alkalinity + 2*precipitate` for fixed-inventory reaction steps, with Na and Cl
reported as reaction-inert inventories. Reaction subcycling can limit accepted
substeps by `maxVmChangePerStep` and records requested time, accepted substep
count, maximum accepted `Vm` change, and the active limiter. The `Vm >= 0.6`
blocked-cell mask and topology-change diagnostics are implemented. Topology
changes now call a pluggable flow recomputation hook; the default path records
only a marked blocked-fraction proxy for relative permeability, pressure drop,
and flow rate. Setting `spec.yoonFlowSolver =
'finite_volume_darcy_pressure'` runs a topology-aware finite-volume Darcy
pressure solve around substrate and blocked cells, then records pressure,
velocity, flux, and provenance fields (`flowSolver`, `flowIsProxy`,
`flowIsStokes`). This is a spatial flow-redistribution diagnostic, not a Stokes
solver. When finite-volume transport is active, these velocity fields now drive
the advective CFL and upwind conservative-component advection; the default
blocked-fraction proxy has no velocity field and does not alter transport. A
lightweight `spec.yoonFlowSolver = 'finite_difference_stokes_brinkman'` path
now solves a cell-centered Stokes-Brinkman saddle-point system on the Yoon grid
and records `flowIsStokes = true`, velocities, pressure, fluxes, and linear
residual diagnostics. This is the first independent Stokes-family hook for the
Yoon micro-continuum path, but it still needs production-scale validation,
grid convergence, and comparison against the final Stokes backend before it can
satisfy the full Yoon flow-feedback benchmark. It does not yet include full
production transport-reaction coupling.

Short diffusion-feedback sensitivity smokes are implemented for `n = 0`, `2`,
and `3`. They are intended to verify that `D_eff = D * (1 - Vm)^n` changes
subsequent mixing and precipitate inventory; they are not a substitute for the
planned 10/5/2.5 um grid-convergence study. The sensitivity runner can write a
`diffusion_feedback_manifest.json` record with the Case 2/1/3 exponent order,
final area, precipitate inventory, mean effective diffusivity, trend checks,
and `diffusionFeedbackAccepted`; readiness audits require that acceptance flag
before quantitative diffusion-feedback claims. The gate also verifies the
underlying exponent sequence, finite nonnegative three-case area, precipitate
inventory, and mean effective diffusivity arrays, all `feedbackCoupled` entries,
and the area/mass trends, so a standalone `diffusionFeedbackAccepted` flag
cannot bypass failed or missing numeric evidence.
`precip_RunYoonGridConvergenceSmoke` now runs the fixed-geometry diagnostic
path on multiple grid sizes and records final-area differences against the
finest smoke grid. It can derive cases from `targetGridSpacing_um`, including
the planned `[10, 5, 2.5]` um sequence, and writes the requested target spacing
to both summary and manifest outputs. The manifest records requested and actual
spacing evidence, including `actualDx_um`, `actualDy_um`, `actualNumX`,
`actualNumY`, `requiredTargetGridSpacing_um`,
`targetGridSpacingSequenceComplete`, `actualGridSpacingSequenceComplete`,
`caseManifestPaths`, `caseManifestFilesVerified`, `gridConvergenceTolerance`,
`gridConvergenceWithinTolerance`, `productionGridConvergenceValidated`,
`gridConvergenceAcceptanceCriteria`, and `gridConvergenceAccepted`. Default
smoke runs remain non-quantitative because
they leave `productionGridConvergenceValidated = false`; this is a
reproducibility scaffold for later convergence work, not evidence that the
final production convergence study has been completed. Readiness audits use
`gridConvergenceAccepted` plus the underlying sequence, tolerance, finite
relative-difference, finite tolerance, and production-validation fields for the
grid gate, so legacy or hand-authored `isQuantitativeBenchmark` flags cannot
bypass the convergence criteria.

Short Case 5 dissolution smokes apply `dissolutionFactor = 300` only after a
center-band blocked cell activates the Case 5 boost, and only to negative
Vaterite rates. This verifies the delayed switch behavior, but it is not by
itself the full 13/18/118 min Case 1 versus Case 5 reproduction.
`precip_RunYoonCase1Case5FixedGeometryComparison` now runs the fixed-geometry
diagnostic path for both Case 1 (`dissolutionFactor = 1`) and Case 5
(`dissolutionFactor = 300`), preserving separate case manifests and a combined
summary CSV/JSON. The comparison manifest now records explicit acceptance
ingredients: target-time completeness, reaction-mass ledger acceptance,
delayed/dissolution-only Case 5 factor provenance, `case5BoostActivated`,
`case5ActivationTime_s`, finite nonnegative Case 1/Case 5 final precipitate
moles and total area, Case 5 final precipitate/area trend checks,
`productionComparisonValidated`, and the derived
`case1Case5Accepted`. This is the planned Case 1/Case 5 output scaffold, still
marked non-quantitative until production transport-reaction coupling, exact
geometry/masks, and grid convergence are complete.

Target-time `Vm` snapshot capture and export are implemented for the
fixed-geometry Yoon smoke path. Snapshot `.mat` files and area CSVs are
diagnostic artifacts until the planned coupled benchmark run and grid
convergence are available.
Static substrate geometry can now be loaded from `spec.substrateMaskFile` as a
`.mat` file with `substrateMask` or as a CSV table with `row`, `col`, and
`substrateMask`. The solver state and fixed-geometry manifest record the
geometry source and file path. This is the ingestion point for future digitized
Zhang/Yoon geometry; the repository still needs the actual calibrated image,
processing script, and versioned mask file.
`precip_LoadZhangYoonGeometryPackage` now validates a complete calibrated
geometry package with source image, calibration CSV, processing script,
substrate mask, region mask file, and uncertainty CSV. Setting
`spec.geometryPackageDir` on fixed-geometry runs attaches the package and
records quantitative-geometry provenance in the manifest. The default
approximate cylindrical-post geometry remains non-quantitative.
Reference curves can now be loaded with `precip_LoadReferenceCurves`, which
validates the required CSV schema, allowed source/case/region combinations, and
finite time/area values. The loader records provenance metadata and marks the
asset as non-quantitative because the current CSV is still approximate visual
digitization rather than a fully archived WebPlotDigitizer/source-data package.
`precip_LoadReferenceDigitizationPackage` now validates that a quantitative
package includes the source screenshot, WebPlotDigitizer project, calibration
CSV, raw export CSV, conversion script, uncertainty CSV, and converted
reference CSV. `precip_LoadReferenceCurves` can attach that package provenance
and only marks complete packages as quantitative.
Yoon area region masks are now configurable through `spec.regionMasks`; the
default remains an x-window fallback. Quantitative comparisons still require
digitized first-pore and first-three-pore masks from the published geometry,
with mask provenance recorded in the manifest. The same masks can also be
loaded from `spec.regionMaskFile` using either `.mat` variables
`firstPoreMask` / `firstThreePoresMask` or a CSV table with `row`, `col`,
`firstPoreMask`, and `firstThreePoresMask` columns.

`precip_RunYoonFixedGeometryBenchmark` provides a manifest-writing diagnostic
runner for the fixed-geometry path, including default 13/18/118 min targets.
Its manifest explicitly marks the output as non-quantitative benchmark
evidence, records the active region-mask source, and stores a center-band
morphology diagnostic from the final active precipitate centroid relative to
the split inlet. It records `capturedTargetTimes_s`, `missingTargetTimes_s`,
`numSnapshots`, and `targetSnapshotsComplete` so target-time evidence is based
on actual snapshot capture rather than requested times alone. The readiness
gate requires all 13/18/118 min targets in `capturedTargetTimes_s`, empty
`missingTargetTimes_s`, and at least three exported snapshots. It also records
flow-feedback audit fields:
`flowTopologyChangedAny`, `totalFlowTopologyChangedSteps`,
`flowRecomputedAfterTopologyChange`, `flowBackendProductionValidated`, and
`flowFeedbackAccepted`, plus finite `finalRelativePermeability`,
`finalPressureDropRelative`, and `finalFlowRateRelative` metrics from the final
recomputed flow field. The fixed-geometry runner derives
`flowBackendProductionValidated` from the final non-proxy Stokes solver
provenance, recorded linear-residual threshold, divergence-residual threshold,
and inlet/outlet boundary-flux closure threshold instead of accepting a manual
option. A `hyphm_stokes` bridge maps Yoon blocked/substrate masks to a
positive-fluid/negative-solid level set, calls an injected HyPHM Stokes solver,
and maps returned velocity samples back to the Yoon transport grid; it errors
if no solver function is supplied. Independent production-flow validation now
writes `yoon_production_flow_validation_manifest.json` for `empty_channel`,
`single_obstacle`, and `blocked_column` cases, and the readiness gate requires
that manifest to pass before accepting the `production_stokes` requirement.
`precip_AuditYoonBenchmarkReadiness` is the conservative quantitative-claim
gate. It requires 13/18/118 min target coverage, complete quantitative geometry
with exported snapshots, complete quantitative geometry and reference packages,
quantitative 10/5/2.5 um grid convergence,
production-validated Stokes flow provenance, accepted reaction mass ledgers,
accepted center-band morphology, accepted `Vm >= 0.6` flow-feedback
recomputation evidence, literature-locked Vaterite rate constants, accepted
PHREEQC/Yoon speciation cross-validation with finite pH/SI differences within
recorded thresholds, accepted no-reaction split-inlet
conservative transport, accepted Case 1/Case 5 comparison evidence, accepted
diffusion-feedback Case 2/1/3 sensitivity evidence, explicit no-finite-clipping
evidence, and a quantitative benchmark manifest carrying required output
evidence paths. Current smoke manifests are expected to fail this audit.
The reference-package gate consumes digitization provenance directly: source
figure, screenshot, WebPlotDigitizer project, calibration, raw export,
conversion script, uncertainty record, converted CSV, package note, and an
empty `missingAssets` record must be present. The gate also requires
`assetFilesVerified = true` from the loader, verified calibration/raw-export
/uncertainty tables, and `numReferenceRows > 0`, so nonempty provenance
strings alone cannot stand in for checked assets, machine-readable
digitization tables, and converted reference rows.
The geometry-package gate consumes fixed-geometry provenance fields directly:
package directory, package name, package note, substrate mask file, and region
mask file must be nonempty, and the note must not indicate missing or
incomplete or invalid geometry assets. It also requires
`geometryPackageAssetFilesVerified = true`, verified geometry calibration and
uncertainty tables, a positive substrate-mask cell count, a positive first-pore
mask cell count, and a first-three-pores mask cell count no smaller than the
first-pore count.
The quantitative-manifest gate consumes output provenance directly: `outputRoot`,
`snapshotDir`, `areaCsv`, `flowDiagnosticsCsv`, `transportDiagnosticsCsv`,
`reactionMassLedgerCsv`, `reactionDiagnosticsCsv`, and a nonempty `matFiles`
snapshot list must accompany `isQuantitativeBenchmark = true`. The gate also
requires `outputEvidenceFilesVerified = true` so listed output evidence has
been checked after export. It additionally requires
`transportMode = finite_volume` and `reactionSubcycling = true`, so analytic
transport smokes cannot be relabeled as quantitative benchmark runs.
The rate-constant gate requires a non-smoke/non-pending source plus finite
nonnegative `yoonRateK1`, `yoonRateK2`, and `yoonRateK3` values with
`yoonRateUnits = mol_cm-2_s-1`, nonempty source DOI/equation fields,
`yoonRateMineralPhase = Vaterite`, and `yoonRateSourceValuesVerified = true`;
the accepted flag alone is not enough.
The flow-feedback gate consumes the fixed-geometry criteria fields directly:
topology changes must be observed, flow recomputation count must cover the
topology-change count, the final flow must be Stokes rather than a proxy, and
the backend must carry production-validation provenance with a nonempty
`finalFlowSolver` name that identifies a Stokes backend rather than Darcy or a
proxy. It also requires finite `finalFlowLinearResidualRelative` no greater
than finite `maxAcceptedFlowLinearResidualRelative`, so solver labels cannot
replace numerical convergence evidence.
The reaction-mass gate likewise consumes the fixed-geometry ledger fields:
finite Ca, C, and alkalinity relative closure errors must be no greater than
`1e-10`, not merely accompanied by a standalone accepted flag.
The center-band morphology gate checks that active precipitate area is positive
and that finite centroid-distance/tolerance fields place the final precipitate
centroid within the accepted distance from the split inlet.
The passive-transport gate consumes the manifest criteria fields directly,
including limiter provenance, minimum component concentration, inert Na/Cl mass
closure, mixing symmetry, rejected substeps, and the recorded acceptance
thresholds. Finite-volume passive transport evidence should come from the
micro-continuum state path, where `componentMoles` is updated by mol/s face
fluxes and boundary closure is checked against cell inventory rather than
plain concentration sums. The accepted evidence must record
`transportMode = finite_volume`, `initialComponentSource = initial`, and finite
`maxAbsInertBoundaryClosureRelativeError` no greater than
`maxAcceptedBoundaryClosureRelativeError`; the inert initial-to-final inventory
change is diagnostic for the open split inlet, not the conservation criterion.
It also requires the manifest `componentNames` to include
`Ca_total`, `C_total`, `Na_total`, `Cl_total`, and `Alkalinity`, with no
`H_total` or free-H transport field; a standalone `passiveTransportAccepted`
flag is not sufficient without those fields.
The Case 1/Case 5 gate consumes the comparison manifest criteria fields
directly, including target-time completeness, reaction-mass acceptance,
dissolution-only and delayed-activation Case 5 factor provenance, finite
nonnegative Case 1/Case 5 final precipitate moles and total area, final Case 5 precipitate/area trend
checks, and `productionComparisonValidated`; a standalone
`case1Case5Accepted` flag is not sufficient without those fields.
`precip_RunYoonBenchmarkReadinessAudit` writes the same readiness result to a
requirements CSV and JSON manifest from fixed-geometry/grid/reference evidence,
plus optional chemistry, passive-transport, Case 1/Case 5, and
diffusion-feedback manifests, so benchmark outputs can carry an explicit
machine-readable pass/fail record.

A transport subcycle controller is available for conservative component
transport candidates. It rejects negative Ca/C/Na/Cl or mass-balance drift and
shrinks the substep rather than clipping concentrations. The passive
split-inlet diagnostic now defaults to finite-volume micro-continuum transport
from DI-water initial components, and the legacy analytic field remains only an
explicit compatibility option. The finite-volume path records a boundary-flux
ledger for open split-inlet runs and a
closed-boundary mass-conservation check for operator testing. The short Case 1
driver can now call this candidate before each reaction step and export
transport diagnostics. The finite-volume driver path covers each requested
macro step with accepted shrink/retry substeps and reports total advanced time,
accepted substep count, rejected-step count, boundary-closure error, stable
transport step, and active advective/diffusive CFL limiter. This remains
smoke-level coupling: it still needs production-scale verification, grid
convergence, exact geometry masks, and real flow recomputation before the
benchmark can be treated as coupled validation.

## Pending Source Lockdown

Before quantitative validation, add the digitized geometry, figure
digitization projects, raw reference CSV exports, unit-conversion scripts, and
digitization uncertainty records used for first-pore and first-three-pore area
comparisons. Geometry source images, calibration files, processing scripts,
mask files, and geometry uncertainty records must be archived as validated
geometry packages.

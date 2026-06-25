# RTM CaCO3 Precipitation Module

This folder contains the precipitation-specific RTM/PNM workflow for the
Zhang 2010 / Yoon 2012 calcium-carbonate micromodel benchmark. It is intended
to keep precipitation experiments isolated from the original dissolution solver
in `ReactiveTransport/RTM/PNM_beauty3.m`.

## Purpose

- Run the Yoon seeded micro-continuum benchmark scaffold with Vaterite as the
  mineral phase and `Vm` as the precipitated-solid state variable.
- Keep the public Zhang/Yoon entrypoint on
  `precip_ZhangYoonBenchmarkSpec` plus
  `precip_RunYoonFixedGeometryBenchmark`.
- Preserve the older signed Calcite level-set solver only as a legacy
  diagnostic path.
- Export local diagnostic area-time curves and 13, 18, and 118 min snapshots
  for comparison with Zhang 2010 experiment 3 and Yoon 2012 Case 1 / Case 5.

## Model Layering

The current benchmark path is deliberately limited to
`modelFamily = yoon_seeded_microcontinuum` and
`precipitationMode = yoon_seeded_microcontinuum`. It represents the Yoon
seeded-growth reproduction in which top/bottom surfaces and micromodel walls
provide the initial reactive area.

Later Deng-style modes such as `deng_homogeneous_nucleation` and
`surface_growth` are deferred until the Yoon quantitative gate is complete.
`precip_YoonMicrocontinuumSolver` fails closed if either mode is requested, so
random nuclei fields such as `nucleiNumberDensity`, `crystalMoles`,
`crystalSurfaceArea`, `nucleationEventCount`, and `randomSeed` do not silently
enter the Zhang/Yoon benchmark state.

## Main Entry Points

Run the local Zhang/Yoon benchmark through the authoritative Yoon/Vm
entrypoint:

```matlab
addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate');
result = run_zhang_yoon_caco3_precipitation_benchmark();
disp(result.manifestPath);
```

Run the legacy signed Calcite level-set diagnostic:

```matlab
legacy = run_legacy_signed_calcite_levelset_diagnostic();
disp(legacy.resultsDir);
```

Run the new Yoon micro-continuum foundation smokes:

```matlab
addpath(genpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate'),'-begin');
spec = precip_ZhangYoonBenchmarkSpec();
state = precip_YoonMicrocontinuumSolver('initialize', spec);
passive = precip_RunPassiveSplitInletBenchmark(spec);
mix = precip_BuildYoonMixingSeries(spec, [1; 0.5; 0]);
chem = precip_YoonCarbonateEquilibrium(mix, spec);
area = precip_ComputeYoonReactiveArea(state, spec);
rate = precip_YoonVateriteRate(chem, spec);
deff = precip_UpdateEffectiveDiffusivity(state, spec);
next = precip_AdvanceReaction(state, spec, 60);
short = precip_RunYoonCase1Short(spec, struct('endTime_s', 13*60, 'dt_s', 30));
```

Run a fixed-geometry diagnostic Yoon benchmark and export target-time `Vm`
snapshots:

```matlab
spec = precip_ZhangYoonBenchmarkSpec();
diagnostic = precip_RunYoonFixedGeometryBenchmark(spec);
disp(diagnostic.manifestPath);
```

Run a fixed-geometry diagnostic Case 1 / Case 5 comparison:

```matlab
spec = precip_ZhangYoonBenchmarkSpec();
comparison = precip_RunYoonCase1Case5FixedGeometryComparison(spec);
disp(comparison.manifestPath);
```

Run a small multi-grid diagnostic smoke and export a convergence summary:

```matlab
spec = precip_ZhangYoonBenchmarkSpec();
convergence = precip_RunYoonGridConvergenceSmoke(spec);
disp(convergence.manifestPath);
```

Build convergence cases from target physical grid spacing, including the
planned 10/5/2.5 um sequence:

```matlab
spec = precip_ZhangYoonBenchmarkSpec();
convergence = precip_RunYoonGridConvergenceSmoke(spec, struct( ...
    'targetGridSpacing_um', [10, 5, 2.5]));
disp(convergence.summary(:, {'label', 'targetGridSpacing_um', 'numX', 'numY'}));
```

Write a quantitative-readiness audit report from benchmark evidence:

```matlab
auditReport = precip_RunYoonBenchmarkReadinessAudit(struct( ...
    'fixedGeometryManifestPath', diagnostic.manifestPath, ...
    'gridConvergenceManifestPath', convergence.manifestPath, ...
    'referenceCsv', fullfile(pwd, 'reference_data', ...
        'zhang_yoon_reference_curves.csv'), ...
    'outputRoot', pwd));
disp(auditReport.audit.requirements);
```

Build only the legacy level-set diagnostic configuration:

```matlab
cfg = precip_ConfigureZhangYoonBenchmark();
```

Run helper and integration tests:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate','-begin'); addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\couplePhreeqc'); results = runtests('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate\tests\test_precip_signed_helpers.m'); assertSuccess(results);"
```

## Key Files

- `run_zhang_yoon_caco3_precipitation_benchmark.m`: user-facing Zhang/Yoon
  Yoon/Vm runner; it calls `precip_ZhangYoonBenchmarkSpec` and
  `precip_RunYoonFixedGeometryBenchmark`.
- `run_legacy_signed_calcite_levelset_diagnostic.m`: legacy signed Calcite
  level-set diagnostic wrapper around `precip_PNM_beauty3`.
- `precip_ConfigureZhangYoonBenchmark.m`: geometry, flow, chemistry, PHREEQC,
  comparison-time, and output configuration for the legacy diagnostic path.
- `precip_PNM_beauty3.m`: copied legacy precipitation solver. It is retained
  for signed Calcite level-set diagnostics, not as the authoritative Yoon/Vm
  benchmark path.
- `precip_RunPhreeqcCalciteBatchSigned.m`: signed PHREEQC batch runner.
- `precip_WritePhreeqcSpeciesTable.m`: local PHREEQC species CSV export,
  including signed calcite amount/rate diagnostics.
- `precip_ExportPhreeqcSpeciesPlots.m`: local PHREEQC concentration and signed
  calcite diagnostic plot export.
- `precip_PrepareConcentrationFaceData.m`: local P0 triangle mask helper used
  by PHREEQC diagnostic plots.
- `precip_ComputeBoundaryInletMassFlux.m`: diagnostic inlet solute flux helper
  for uniform and split-left inlet cases.
- `precip_ComputeBoundaryInletConcentrationScale.m`: diagnostic inlet
  concentration-scale helper for stability flags.
- `precip_ComputeSignedCalciteInterfaceRatePerArea.m`: converts signed calcite
  delta to interface rate per area.
- `precip_CreateTransportMultiInlet.m`: split inlet boundary helper.
- `precip_ComputePrecipitationAreaMetrics.m`: local area metrics.
- `precip_ExportBenchmarkSnapshots.m`: stable snapshot export for comparison
  times.
- `precip_CompareZhangYoonBenchmark.m`: comparison CSV/PNG/report generator.
- `reference_data/zhang_yoon_reference_curves.csv`: source-bounded approximate
  digitized reference rows.
- `benchmark/precip_ZhangYoonBenchmarkSpec.m`: source-bounded constants for
  the Vaterite/Vm Yoon micro-continuum path.
- `benchmark/precip_LoadZhangYoonGeometry.m`: loads static substrate masks
  from `.mat` or row/column `.csv` files so future digitized Zhang/Yoon
  geometry can replace the approximate post array without changing `Vm`
  semantics.
- `benchmark/precip_LoadZhangYoonGeometryPackage.m`: validates a calibrated
  geometry/mask package containing the source micromodel image, calibration
  CSV, mask-processing script, substrate mask, first-pore/first-three-pore
  region masks, and uncertainty CSV. A complete package can be attached with
  `spec.geometryPackageDir`.
- `benchmark/precip_LoadReferenceCurves.m`: provenance-checked loader for
  repository or user-supplied Zhang/Yoon digitized reference curves. It
  validates source/case/region combinations, rejects header-only files, fills
  normalized values from Yoon cm2 rows by group, and marks the asset as
  non-quantitative until source-lockdown files are added.
- `benchmark/precip_LoadReferenceDigitizationPackage.m`: validates a
  source-locked digitization package containing the source screenshot,
  WebPlotDigitizer project, axis calibration CSV, raw export CSV, conversion
  script, uncertainty CSV, and converted reference CSV. `precip_LoadReferenceCurves`
  can use this package provenance to mark complete packages as quantitative.
- `solver/precip_YoonMicrocontinuumSolver.m`: initializes static substrate,
  dynamic `Vm`, blocked-mask, and conservative component state without
  level-set or Calcite inventory assumptions.
- `solver/precip_RunPassiveSplitInletBenchmark.m`: no-reaction split-inlet
  conservative transport smoke for Ca/C/Na/Cl/alkalinity fields. The default
  path starts from DI-water initial components and advances real
  micro-continuum finite-volume split-inlet steps with reaction disabled. It
  can write a `passive_transport_manifest.json` acceptance record with
  nonnegativity, inert boundary-closure, mixing-symmetry, and substep
  diagnostics. It records the real initial-to-final inventory change separately
  from the split-inlet boundary-flux closure ledger.
- `solver/precip_AdvanceConservativeTransport2D.m`: explicit no-clipping
  finite-volume smoke transport step for conservative components, including
  closed-boundary mass conservation checks and split-inlet boundary-flux
  ledger diagnostics. It supports substrate masks, blocked-cell advection
  masks, and local `effectiveDiffusivity_cm2_s` face fluxes.
- `solver/precip_AdvanceMicrocontinuumTransport2D.m`: state-based
  finite-volume transport step whose conservative variable is
  `state.componentMoles`. It assembles mol/s face fluxes, updates inventories,
  and then recovers water-phase `components`/`aqueousConcentration` from
  `fluidVolumeFraction * cellVolume`.
- `solver/precip_AssembleComponentFaceFluxes.m`: assembles component face
  fluxes for micro-continuum transport, with static substrate cells blocking
  advection and diffusion, blocked `Vm` cells blocking advection, and local
  `effectiveDiffusivity_cm2_s` controlling diffusive faces.
- `solver/precip_ComputeMicrocontinuumStableDt.m`: state-aware CFL wrapper
  that passes substrate, blocked, velocity, and local diffusivity fields into
  the transport stable-step calculation.
- `solver/precip_RefreshYoonComponentMolesFromAqueous.m` and
  `solver/precip_RefreshYoonAqueousFromComponentMoles.m`: synchronize the
  Yoon state between conservative cell inventories and water-phase
  concentrations. `componentMoles` is the mass-audit state; concentrations are
  used for chemistry and output.
- `solver/precip_ComputeTransportStableDt.m`: computes explicit transport
  stable step limits from advective and diffusive CFL controls for the
  finite-volume Yoon path, including local maximum effective diffusivity when
  supplied by the Yoon state.
- `solver/precip_RunTransportSubcycle.m`: reject/shrink/retry controller for
  conservative transport candidates; rejects negative components and
  mass-balance drift instead of clipping concentrations.
- `solver/precip_AdvanceReaction.m`: fixed-geometry reaction step chaining
  carbonate speciation, Vaterite rate evaluation, reactive area, `Vm`, and
  effective diffusivity updates.
- `solver/precip_ComputeYoonReactionStableDt.m`: computes a reaction substep
  limit from `maxVmChangePerStep` and the instantaneous `Vm` rate.
- `solver/precip_RunYoonReactionSubcycle.m`: advances Yoon reaction over a
  macro step using `Vm`-limited reaction substeps and records subcycle
  diagnostics.
- `solver/precip_RunYoonCase1Short.m`: fixed-geometry, analytic split-inlet
  Yoon Case 1 smoke driver with area and mass diagnostics. It can also run
  an optional `transportMode = 'finite_volume'` path that advances
  conservative components before each reaction step and records transport
  subcycle diagnostics. When the active flow field exposes velocity arrays,
  finite-volume transport uses those velocities for CFL and upwind advection.
  It records whether `reactionSubcycling` was used.
- `solver/precip_UpdateFlowMask.m`: updates Yoon blocked-cell masks using
  `Vm >= 0.6` and reports topology-change diagnostics.
- `solver/precip_RecomputeYoonFlowField.m`: flow recomputation adapter called
  after blocked-mask topology changes. Without a supplied `flowSolverFcn` or
  `spec.yoonFlowSolver` override, it returns a marked proxy diagnostic from
  blocked fraction rather than solving flow.
- `solver/precip_SolveYoonDarcyFlow2D.m`: topology-aware finite-volume Darcy
  pressure recomputation for `spec.yoonFlowSolver =
  'finite_volume_darcy_pressure'`. It outputs pressure, x/y velocity
  diagnostics, fluxes, relative permeability, and explicit `isStokes = false`
  provenance; final benchmark claims should use the validated Stokes-family
  path instead of this Darcy diagnostic.
- `solver/precip_SolveYoonStokesBrinkman2D.m`: lightweight cell-centered
  finite-difference Stokes-Brinkman recomputation for `spec.yoonFlowSolver =
  'finite_difference_stokes_brinkman'`. It solves a Yoon-grid saddle-point
  system with no-flow substrate/blocked cells and records velocity, pressure,
  flux, residual, and explicit `isStokes = true` provenance. It is an
  auditable Stokes-family diagnostic hook, not yet a production-validated
  replacement for the final benchmark Stokes backend.
- `solver/precip_MapCellMaskToLevelSet.m`: converts Yoon `substrateMask`,
  `blockedMask`, and `Vm >= blockedVmThreshold` cells into a
  positive-fluid/negative-solid level-set package for HyPHM-side Stokes calls.
- `solver/precip_MapHyPHMFluxToCartesianFaces.m`: maps injected HyPHM velocity
  samples back to Yoon cell-centered `velocityX_cm_s`/`velocityY_cm_s` fields
  and records inlet/outlet flux closure plus divergence diagnostics.
- `solver/precip_SolveYoonHyPHMStokes.m`: fail-closed bridge for
  `spec.yoonFlowSolver = 'hyphm_stokes'`. It requires an injected
  `spec.hyphmStokesSolverFcn` or `options.hyphmStokesSolverFcn`, maps the Yoon
  mask to a level set, calls that solver, and maps returned velocity samples
  back to the Yoon finite-volume transport grid.
- `solver/precip_RunYoonProductionFlowValidation.m`: runs independent
  `hyphm_stokes` validation cases (`empty_channel`, `single_obstacle`,
  `blocked_column`) and writes a production-flow manifest with linear
  residual, divergence residual, inlet/outlet closure, and permeability bounds.
- `solver/precip_RunYoonDiffusionFeedbackSensitivity.m`: runs short
  fixed-geometry Case 2/1/3 smokes for `n = 0/2/3` and reports area, total
  precipitate inventory, and mean `D_eff`. It can also write a
  `diffusion_feedback_manifest.json` acceptance record for readiness audits.
- `solver/precip_RunYoonCase5ShortComparison.m`: compares short Case 1 and
  Case 5 dissolution from the same precipitated state, applying
  `dissolutionFactor = 300` only to undersaturated dissolution rates after the
  center-band blocked-cell activation condition is met.
- `solver/precip_RunYoonFixedGeometryBenchmark.m`: fixed-geometry diagnostic
  runner that captures default 13/18/118 min `Vm` snapshots, exports `.mat`
  files and area metrics, and writes `yoon_fixed_geometry_manifest.json`. The
  manifest records a center-band morphology diagnostic from the final active
  precipitate centroid relative to the split inlet. It also records
  flow-feedback evidence fields for `Vm >= 0.6` topology changes, flow
  recomputation, production-validation provenance, and `flowFeedbackAccepted`.
  Geometry package evidence includes `geometryCalibrationTableVerified` and
  `geometryUncertaintyTableVerified`, so placeholder calibration files cannot
  satisfy the quantitative geometry gate.
  Target-time evidence includes `capturedTargetTimes_s`,
  `missingTargetTimes_s`, `numSnapshots`, and `targetSnapshotsComplete`, so
  readiness can distinguish requested 13/18/118 min targets from actually
  exported snapshots. If `isQuantitativeBenchmark = true` is requested, the
  runner requires `transportMode = 'finite_volume'` and
  `reactionSubcycling = true`; explicit analytic transport or disabled
  reaction subcycling is rejected before the run starts.
- `solver/precip_RunYoonCase1Case5FixedGeometryComparison.m`: runs the
  fixed-geometry diagnostic path for Case 1 (`dissolutionFactor = 1`) and Case
  5 (`dissolutionFactor = 300`) into separate case folders, then writes
  `yoon_case1_case5_comparison_summary.csv` and
  `yoon_case1_case5_comparison_manifest.json`. The manifest records
  acceptance-criteria fields for target-time coverage, reaction-mass ledgers,
  delayed/dissolution-only Case 5 factor provenance, activation time, final
  Case 5 precipitate/area trends, and production comparison validation.
- `solver/precip_RunYoonGridConvergenceSmoke.m`: multi-grid diagnostic wrapper
  around the fixed-geometry runner. It writes
  `yoon_grid_convergence_summary.csv` and
  `yoon_grid_convergence_manifest.json`. It accepts explicit `gridCases` or
  derives cases from `targetGridSpacing_um`, including the planned
  `[10, 5, 2.5]` um sequence. The manifest records target-spacing sequence
  completeness, convergence tolerance, production-validation provenance, and
  `gridConvergenceAccepted`.
- `diagnostics/precip_ExportYoonSnapshots.m`: exports captured `Vm` snapshots
  as `.mat` files plus a Yoon area-metrics CSV.
- `diagnostics/precip_BuildYoonRegionMasks.m`: builds first-pore and
  first-three-pore masks for Yoon area metrics, using explicit `spec.regionMasks`
  when provided and a documented x-window fallback otherwise.
- `diagnostics/precip_LoadYoonRegionMasks.m`: loads first-pore and
  first-three-pore masks from `.mat` files or CSV row/column tables for
  future digitized-geometry inputs.
- `chemistry/precip_YoonCarbonateEquilibrium.m`: zero-dimensional carbonate
  equilibrium smoke used to verify that pH 10.9 is preserved and Vaterite
  supersaturation appears in mixed cells.
- `chemistry/precip_SpeciateYoonComponents.m`: chemistry backend dispatcher
  for the Yoon path. The default `yoon_equilibrium` backend provides aqueous
  speciation; `iphreeqc_speciation` is an explicit injection interface via
  `spec.iphreeqcSpeciationFcn` and still fails clearly when no PHREEQC-backed
  function is supplied.
- `chemistry/precip_IPhreeqcSpeciation.m`: zero-dimensional PHREEQC aqueous
  speciation backend for Yoon components. It builds `SOLUTION` blocks from
  `Ca_total`, `C_total`, `Na_total`, `Cl_total`, and `Alkalinity`, runs a
  persistent `rtm.phreeqc.PhreeqcSession`, and parses pH, activities, and
  Vaterite SI. The input defines a `PHASES` entry for Vaterite from
  `spec.vateriteKsp` because the bundled `phreeqc.dat` does not include that
  phase and preserves signed Alkalinity values instead of clipping acidic
  samples to zero. Unit tests use a mock engine; a local three-point IPhreeqcCOM
  smoke has run, but this is not yet full quantitative PHREEQC/Yoon acceptance.
- `chemistry/precip_RunYoonSpeciationMixingSeries.m`: builds the
  `f = 0:0.01:1` zero-dimensional mixing series, compares local Yoon carbonate
  equilibrium against an optional PHREEQC speciation backend, and reports pH/SI
  differences plus peak supersaturation locations. It can also write a
  `yoon_speciation_mixing_manifest.json` acceptance record for readiness
  audits, including 101-point mixing-fraction coverage evidence.
- `chemistry/precip_YoonVateriteRate.m`: signed Yoon/Chou-style Vaterite rate
  evaluator. Defaults are locked to the Yoon/Chou Vaterite constants in
  `mol_cm-2_s-1`.
- `precipitation/precip_ComputeYoonReactiveArea.m`: Yoon reactive area for
  fluid cells, including top/bottom area, vertical faces adjacent to static
  substrate, and exposed vertical faces adjacent to fully filled `Vm = 1`
  precipitate cells.
- `precipitation/precip_UpdateVateriteVolumeFraction.m`: stoichiometric
  Ca/C/alkalinity and dynamic `Vm` update. Precipitation is limited by Ca/C and
  available `Vm`, not by a nonnegative Alkalinity assumption. The update
  subtracts accepted Vaterite moles from `state.componentMoles` and then
  recomputes water-phase concentration using the reduced fluid volume.
- `precipitation/precip_UpdateEffectiveDiffusivity.m`: computes
  `D_eff = D * (1 - Vm)^n`.
- `diagnostics/precip_ComputeYoonAreaMetrics.m`: Yoon-style projected
  precipitated area using `Vm > 0.05` and optional region masks.
- `diagnostics/precip_ComputeYoonReactionMassLedger.m`: Yoon reaction
  stoichiometry ledger for `Ca_total + precipitate`, `C_total + precipitate`,
  `Alkalinity + 2*precipitate`, and inert Na/Cl inventories.
- `diagnostics/precip_AuditYoonBenchmarkReadiness.m`: conservative
  quantitative-claim gate for the Zhang/Yoon path. It checks 13/18/118 min
  target coverage, quantitative geometry package, quantitative reference
  package, production 10/5/2.5 um grid-convergence evidence, production
  Stokes flow provenance, reaction mass acceptance, center-band morphology
  acceptance, accepted `Vm >= 0.6` flow-feedback recomputation,
  literature-locked Vaterite rate constants, accepted PHREEQC/Yoon speciation
  cross-validation, accepted passive split-inlet conservative transport,
  accepted Case 1/Case 5 comparison evidence, accepted
  diffusion-feedback Case 2/1/3 sensitivity evidence, explicit
  no-finite-clipping evidence, and quantitative manifest status with required
  output evidence paths.
- `diagnostics/precip_RunYoonBenchmarkReadinessAudit.m`: reads benchmark
  manifest/reference evidence, runs the readiness gate, and writes
  `yoon_benchmark_readiness_requirements.csv` plus
  `yoon_benchmark_ladder.csv` and
  `yoon_benchmark_readiness_manifest.json`. The ladder report maps the
  readiness evidence onto B0-B8 and marks later stages as not entered after
  the first failed stage.

## Output Files

Each run writes to:

```text
outputs/rtm_runs/<runName>/
```

Important files include:

- `run_metadata.json`: configuration and output manifest.
- `global_evolution_log.csv`: porosity, permeability, surface area, and other
  global RTM quantities.
- `stability_diagnostics_log.csv`: concentration range, mass-balance residual,
  local velocity, CFL, reaction CFL, and diagnostic flags.
- `phreeqc_results/phreeqc_summary_log.csv`: PHREEQC pH, saturation index, and
  signed calcite summary statistics.
- `precipitation_area_timeseries.csv`: local precipitation/net-solid area
  metrics.
- `precipitation_area_timeseries.png`: area-time diagnostic plot.
- `benchmark_comparison_times_log.csv`: records captured 13, 18, and 118 min
  comparison states.
- `benchmark_snapshots/benchmark_snapshot_013min.png`,
  `benchmark_snapshot_018min.png`, and `benchmark_snapshot_118min.png`:
  exported interface snapshots when the run reaches those times.
- `zhang_yoon_area_comparison.csv`, `zhang_yoon_area_comparison.png`, and
  `benchmark_comparison_report.md`: comparison outputs generated by
  `precip_CompareZhangYoonBenchmark`.
- `yoon_flow_diagnostics.csv`: diagnostic Yoon-path table written by
  `precip_RunYoonFixedGeometryBenchmark`, including blocked-cell counts,
  flow recomputation counts, relative permeability, pressure-drop proxy, and
  relative flow-rate proxy.
- `yoon_transport_diagnostics.csv`: diagnostic Yoon-path transport table
  written by `precip_RunYoonFixedGeometryBenchmark`, including candidate name,
  requested and accepted transport substeps, total advanced time,
  accepted-substep count, rejected-step count, boundary-closure error, stable
  transport step, and active stability limiter.
- `yoon_reaction_mass_ledger.csv`: diagnostic Yoon-path reaction ledger written
  by `precip_RunYoonFixedGeometryBenchmark`, including relative closure errors
  for total Ca, total C, alkalinity equivalents, and reaction-inert Na/Cl.
- `yoon_reaction_diagnostics.csv`: diagnostic Yoon-path reaction-step table
  written by `precip_RunYoonFixedGeometryBenchmark`, including requested and
  accepted reaction time, total advanced time, accepted substep count, maximum
  accepted `Vm` change, stable reaction step, and active limiter.
- Fixed-geometry manifests also record `yoonRateK1`, `yoonRateK2`,
  `yoonRateK3`, `yoonRateUnits`, `yoonRateSource`, and
  `yoonRateConstantsLocked`, plus `yoonRateSourceDoi`,
  `yoonRateSourceEquation`, `yoonRateMineralPhase`, and
  `yoonRateSourceValuesVerified` so Vaterite rate provenance is
  machine-readable.
- `yoon_grid_convergence_summary.csv` and
  `yoon_grid_convergence_manifest.json`: optional multi-grid smoke outputs
  from `precip_RunYoonGridConvergenceSmoke`, including each case manifest path
  the requested `targetGridSpacing_um` when provided, and the relative
  final-area difference from the finest smoke grid. The manifest also records
  `requiredTargetGridSpacing_um`, `targetGridSpacingSequenceComplete`,
  `gridConvergenceTolerance`, `gridConvergenceWithinTolerance`,
  `productionGridConvergenceValidated`, `gridConvergenceAcceptanceCriteria`,
  and `gridConvergenceAccepted`.
- `yoon_case1_case5_comparison_summary.csv` and
  `yoon_case1_case5_comparison_manifest.json`: optional fixed-geometry Case 1
  / Case 5 comparison outputs from
  `precip_RunYoonCase1Case5FixedGeometryComparison`, including per-case
  manifest paths, final areas, precipitate inventory, flow solver provenance,
  reaction-mass acceptance flags, `case1Case5Accepted`, and dissolution-only
  Case 5 factor provenance. The manifest also records
  `case1Case5AcceptanceCriteria`, `caseTargetTimesComplete`,
  `caseReactionMassAcceptedAll`, `case5FinalPrecipitateMolesLessThanCase1`,
  `case5FinalAreaNoGreaterThanCase1`, `case1FinalPrecipitateMoles`,
  `case5FinalPrecipitateMoles`, `case1FinalTotalPrecipitatedArea_cm2`,
  `case5FinalTotalPrecipitatedArea_cm2`, and
  `productionComparisonValidated`.
- `yoon_speciation_mixing_manifest.json`: optional chemistry cross-validation
  manifest from `precip_RunYoonSpeciationMixingSeries`, recording pH/SI
  differences, acceptance thresholds, PHREEQC backend provenance,
  `numSamples`, `mixingFractionStep`, `mixingFractionsCover101`, and
  `chemistrySpeciationAccepted`.
- `passive_transport_manifest.json`: optional no-reaction split-inlet transport
  manifest from `precip_RunPassiveSplitInletBenchmark`, recording limiter use,
  transport mode, initial component source, minimum component concentration,
  inert inventory change, inert boundary-closure relative error, mixing
  symmetry, rejected substeps, acceptance thresholds, and
  `passiveTransportAccepted`.
- `diffusion_feedback_manifest.json`: optional Case 2/1/3 diffusion-feedback
  sensitivity manifest from `precip_RunYoonDiffusionFeedbackSensitivity`,
  recording the `n = 0/2/3` sequence, final areas, precipitate inventory,
  effective diffusivity, trend checks, and `diffusionFeedbackAccepted`.
- `yoon_production_flow_validation_manifest.json`: optional production Stokes
  validation manifest from `precip_RunYoonProductionFlowValidation`, recording
  required case coverage, HyPHM Stokes provenance, maximum linear/divergence
  residuals, inlet/outlet boundary-closure error, permeability bounds, and
  `productionFlowValidationAccepted`. Readiness audits require this manifest
  for the `production_stokes` gate.
- `yoon_benchmark_readiness_requirements.csv`,
  `yoon_benchmark_ladder.csv`, and
  `yoon_benchmark_readiness_manifest.json`: optional readiness audit outputs
  from `precip_RunYoonBenchmarkReadinessAudit`, listing each quantitative
  benchmark gate, the ordered B0-B8 stage status, and failed requirement or
  stage IDs.

## Signed Calcite Interpretation

PHREEQC selected output is parsed without discarding the sign of calcite
change. The precipitation path uses the following convention:

```text
calcite_delta_moles > 0  -> calcite precipitation, solid grows
calcite_delta_moles < 0  -> calcite dissolution, solid retreats
```

`precip_ComputeSignedCalciteInterfaceRatePerArea.m` converts this signed amount
to a rate per interface area and time step. The copied solver then maps the
rate to level-set normal velocity through the calcite molar volume.

## Current Benchmark Status

The output chain has been exercised through 13, 18, and 118 min with a clipped
diagnostic run:

```text
outputs/rtm_runs/zhang_yoon_caco3_25mM_20260624_095000
```

That run generated the benchmark snapshots, area-time curve, and Zhang/Yoon
comparison report. It should not be treated as quantitative validation because
it used a finite PHREEQC transport limiter:

```text
phreeqcTransportMaxFactor = 2
```

and retained stability flags:

```text
overshoot_c;advective_cfl_gt_1;mass_balance_drift
```

A physically stable validation run still needs to remove persistent CFL and
mass-balance flags without relying on finite clipping as the main stabilizer.

Additional hard-CFL short smokes showed that `allowAdaptiveTimeGridExtension`
can keep appending internal time steps after CFL/concentration shrinkage, and
that H/Ca/C/Na/Cl transport variables are extended together. This fixes the
mechanical "time grid exhausted before endTime" path. It does not solve the
underlying chemistry/transport stability problem: with
`phreeqcTransportMaxFactor = Inf`, the hard-CFL 1 s smoke still developed H
runaway near `t = 0.336 s` and PHREEQC stopped with a water-activity convergence
error. Treat this as evidence that adaptive extension works, not as evidence of
a stable no-clipping precipitation benchmark.

The split-left inlet stability diagnostic now uses the same `inletA` /
`inletB` segment concentrations and `yRange` geometry as the transport
boundary condition when computing the diagnostic inlet solute flux. The
concentration-CFL limiter also uses the split-inlet concentration scale instead
of the legacy premixed scalar. This corrects diagnostic and time-step-control
scales; it does not by itself remove the need for a stable no-clipping
benchmark run.

The signed Calcite PHREEQC builder now preserves positive strong-base inlet
hydrogen concentrations by default (`minHForPHMolL = 1e-12`), so a Zhang
Na2CO3 inlet pH of 10.9 is no longer silently written as pH 7. Nonpositive
legacy H fields still fall back to a neutral pH guard rather than being
interpreted as strong base. This fixes an input-builder truncation bug, but the
signed Calcite/level-set route remains a diagnostic surface-growth path rather
than the Yoon Vaterite micro-continuum benchmark.

The legacy level-set driver also initializes the H+ transport field from
`initialHydrogenConcentration` instead of hard-coding zero. That restores the
configured initial solution for old diagnostic runs; it does not change the
Yoon recommendation to use conservative component transport rather than free
H+ as an independent transported variable.

The Yoon micro-continuum foundation path now exists separately from the
signed Calcite/level-set path. It currently covers benchmark constants,
static-substrate/dynamic-`Vm` initialization, no-reaction split-inlet component
fields, conservative mass-ledger diagnostics, and single-cell carbonate
speciation smokes. It also has unit-tested core reaction operators for
top/bottom plus vertical exposed-face reactive area, signed Vaterite rate sign,
stoichiometric `Vm` updates, effective diffusivity feedback, and a
fixed-geometry reaction step.
The Yoon smoke path stores conservative component fields consistently with the
finite-volume transport and mass ledgers, so the reaction update subtracts
Vaterite stoichiometric amounts from bulk cell component inventories while
`fluidVolumeFraction` records pore-volume feedback. A dedicated reaction mass
ledger checks `Ca_total + precipitate`, `C_total + precipitate`, and
`Alkalinity + 2*precipitate` closure for fixed-inventory reaction steps.
Default Yoon/Chou constants are now source-locked in centimeter units
(`k1 = 8.9e-5`, `k2 = 5.01e-8`, `k3 = 6.6e-11 mol_cm-2_s-1`), with
`Ksp = 1.832e-8 mol2_L-2` and `vateriteMolarVolume_cm3_mol = 37.47`.
Alkalinity is treated as a signed conservative component: PHREEQC input keeps
negative values, and the reaction limiter no longer treats Alkalinity as a
nonnegative elemental inventory.
The reaction step now calls `precip_SpeciateYoonComponents`, making aqueous
speciation an explicit backend choice. `yoon_equilibrium` is the current default;
`iphreeqc_speciation` can call an injected `spec.iphreeqcSpeciationFcn` that
returns PHREEQC-style pH, activities, and Vaterite saturation metadata. The
provided `precip_IPhreeqcSpeciation` implementation uses a persistent
`rtm.phreeqc.PhreeqcSession` with `RunString` and PHREEQC selected output. It
also adds an input-local `PHASES Vaterite` definition from `spec.vateriteKsp`
for reproducibility with the bundled `phreeqc.dat`. It is covered by
mock-session tests, and a local three-point IPhreeqcCOM smoke produced a
supersaturated mixed sample (`SI_Vaterite` about 2.51 with the current
component/pH constraints). PHREEQC-vs-`yoon_equilibrium` quantitative acceptance
over the full mixing series remains pending. A 101-point local IPhreeqcCOM
mixing-series smoke runs through `precip_RunYoonSpeciationMixingSeries`; in the
current setup PHREEQC peaks near `fractionInletA = 0.45` with `SI_Vaterite`
about 2.51, while the local Yoon equilibrium peak is at 0.5. This is a
diagnostic discrepancy to resolve, not a benchmark acceptance. Readiness now
requires the chemistry manifest to record the full `0:0.01:1` fraction grid
through `numSamples = 101`, `mixingFractionStep = 0.01`, and
`mixingFractionsCover101 = true`. The dispatcher still errors clearly when the
function handle is absent.
The Case 1 smoke driver can optionally use `reactionSubcycling = true` to cover
each macro step with substeps limited by `maxVmChangePerStep`; it writes
reaction diagnostics through the fixed-geometry runner. This turns the former
per-cell `Vm` cap into an auditable time-step policy for smoke runs.
The fixed-geometry path is still not a coupled 13/18/118 min Case 1 or Case 5
solver with Stokes flow recomputation.
Short diffusion-feedback sensitivity smokes cover Yoon Case 2 (`n = 0`),
Case 1 (`n = 2`), and Case 3 (`n = 3`) at the operator/driver level. These
are trend checks, not grid-converged benchmark claims. The sensitivity runner
can write `diffusion_feedback_manifest.json`; readiness requires
`diffusionFeedbackAccepted = true` before any quantitative diffusion-feedback
claim is allowed.
The short Case 1 driver now detects `Vm >= 0.6` blocked-cell topology changes
and calls a pluggable `flowSolverFcn(state, spec)` hook when new blockages
appear. It can also use `spec.yoonFlowSolver = 'finite_volume_darcy_pressure'`
to solve a topology-aware finite-volume Darcy pressure problem with left/right
Dirichlet pressure and no-flow top/bottom boundaries, or
`spec.yoonFlowSolver = 'finite_difference_stokes_brinkman'` to solve a small
cell-centered Stokes-Brinkman system on the Yoon grid. The
`spec.yoonFlowSolver = 'hyphm_stokes'` bridge now exists for injected HyPHM
Stokes solvers and errors if no solver function is supplied. Flow diagnostics and the
fixed-geometry manifest now record `flowSolver`, `flowIsProxy`, and
`flowIsStokes`. If no hook or solver override is supplied, the driver records a
clearly marked blocked-fraction proxy. In `transportMode = 'finite_volume'`,
velocity fields from the spatial solvers drive the advective transport CFL and
upwind component advection; the blocked-fraction proxy has no velocity field and
therefore does not alter finite-volume transport. The Stokes-Brinkman path is a
tested diagnostic hook; it still needs production validation and grid
convergence before supporting quantitative Zhang/Yoon claims.
Short Case 5 dissolution smokes verify that `dissolutionFactor = 300` is
applied only when the center band has first blocked and the local Vaterite rate
is negative. Before activation, Case 5 uses the same dissolution rate as Case 1;
after activation, it loses precipitate inventory faster from the same initial
blocked precipitated state.
`precip_RunYoonCase1Case5FixedGeometryComparison` now packages fixed-geometry
Case 1 and Case 5 diagnostic runs side by side, using the same target-time
snapshot/export machinery as the Case 1 runner and recording the per-case
manifest paths plus final area and inventory summaries. Its manifest exposes
the acceptance ingredients rather than only a hard-coded flag: target-time
coverage, reaction-mass acceptance, dissolution-only Case 5 factor provenance,
Case 5 boost activation/time provenance, final precipitate/area trend checks, and
`productionComparisonValidated`. This creates the planned Case 1/Case 5 output
scaffold, but it remains non-quantitative until the coupled production run,
exact geometry, digitized masks, and grid convergence are available.
The short fixed-geometry driver can capture requested target times such as
13, 18, and 118 min and export `Vm` fields plus area metrics. These files are
diagnostic Yoon-path outputs, not final benchmark evidence until the coupled
transport, flow recomputation, and grid-convergence steps are complete.
`precip_RunYoonFixedGeometryBenchmark` packages this diagnostic path into a
single reproducible runner and manifest. It also writes
`yoon_flow_diagnostics.csv` so blocked-mask feedback diagnostics remain
machine-readable alongside `Vm` snapshots. The fixed-geometry manifest now
records which target times were actually captured and whether the requested
target snapshot set is complete; the readiness gate rejects runs that only list
13/18/118 min as requested targets without complete snapshot evidence in
`capturedTargetTimes_s`, empty `missingTargetTimes_s`, and `numSnapshots >= 3`.
Static substrate geometry can now come from `spec.substrateMaskFile`, using a
`.mat` file with `substrateMask` or a CSV with `row`, `col`, and
`substrateMask` columns. The fixed-geometry manifest records
`substrateGeometrySource` and `substrateMaskFile`; the default source remains a
documented approximate cylindrical-post array.
Calibrated geometry packages can now be validated with
`precip_LoadZhangYoonGeometryPackage` and attached to fixed-geometry runs by
setting `spec.geometryPackageDir`. Complete packages supply both the substrate
mask and the first-pore/first-three-pore region masks, and the fixed-geometry
manifest records `geometryPackageDir`, `geometryPackageIsQuantitative`,
`geometryPackageName`, `geometryPackageNote`,
`geometryPackageAssetFilesVerified`, and nonzero substrate/region mask cell
counts.
Reference curves now route through `precip_LoadReferenceCurves`, so the
standalone benchmark asset loader and `precip_CompareZhangYoonBenchmark` share
the same allowed Zhang/Yoon source/case/region schema and missing-normalized
area handling. The repository CSV is still approximate visual digitization.
Source-locked digitization packages can now be validated with
`precip_LoadReferenceDigitizationPackage`; complete packages require the source
screenshot, WebPlotDigitizer project, calibration CSV, raw export, conversion
script, uncertainty CSV, and converted reference CSV before the reference loader
will mark them as quantitative.
First-pore and first-three-pore area diagnostics now route through
`precip_BuildYoonRegionMasks`. Supplying `spec.regionMasks.firstPoreMask` and
`spec.regionMasks.firstThreePoresMask` lets future digitized Zhang/Yoon geometry
masks override the default x-window fallback; mask source metadata is written to
the fixed-geometry manifest. Alternatively set `spec.regionMaskFile` to a `.mat`
file containing `firstPoreMask` and `firstThreePoresMask` variables, or to a
CSV with columns `row`, `col`, `firstPoreMask`, and `firstThreePoresMask`.
The manifest records both `regionMaskSource` and `regionMaskFile`.
The passive split-inlet smoke now exposes transport-subcycle diagnostics. The
default path still uses the analytic split-inlet field for fast smokes, and an
optional `transportMode = 'finite_volume'` path uses the explicit conservative
transport candidate with a split-inlet boundary-flux ledger. This is a
unit-tested transport operator and diagnostic path, not yet the full coupled
transport-reaction benchmark loop.
Finite-volume transport now uses state-derived substrate masks, blocked-cell
advection masks, and local `effectiveDiffusivity_cm2_s` in the actual face
fluxes and CFL estimate; the passive finite-volume ledger compares the initial
and final inventories instead of comparing the final state with itself.
The short Case 1 driver can now run the same `transportMode = 'finite_volume'`
candidate inside its time loop before the reaction step and exports
`yoon_transport_diagnostics.csv` through the fixed-geometry runner. The
finite-volume path covers each requested macro step with accepted shrink/retry
substeps and records total advanced time, accepted substep count, rejected-step
count, boundary-closure error, stable CFL-limited transport step, and active
stability limiter. If a recomputed flow field includes velocity arrays, those
velocities are used for finite-volume advection and the advective CFL. This is
now a driver-level connection of conservative transport to reaction and Darcy
flow diagnostics; it still needs production-scale verification, grid
convergence, exact geometry masks, literature-locked rates, and a
production-validated Stokes backend before quantitative Zhang/Yoon claims.
`precip_RunYoonGridConvergenceSmoke` now packages repeated fixed-geometry runs
on multiple grid sizes and records the final-area spread against the finest
smoke grid. It can now derive cases from physical target spacing, so the
planned 10/5/2.5 um sequence is represented directly in summary and manifest
metadata. This makes convergence audits reproducible and machine-readable, but
the runner only marks convergence accepted when the required 10/5/2.5 um
sequence is present, the finest-grid relative difference is within the recorded
tolerance, and `productionGridConvergenceValidated = true`. Smoke runs leave
that production-validation flag false until the production 10/5/2.5 um study is
actually executed with exact geometry and validated coupling. The readiness
gate consumes both `gridConvergenceAccepted` and the underlying criteria fields:
`targetGridSpacingSequenceComplete`, `actualGridSpacingSequenceComplete`,
`gridConvergenceWithinTolerance`, `productionGridConvergenceValidated`, finite
`maxRelativeTotalAreaDifferenceFromFinest`, and finite
`gridConvergenceTolerance`. It also requires actual exported grid evidence:
`actualDx_um`, `actualDy_um`, `actualNumX`, and `actualNumY` must show the
10/5/2.5 um sequence and increasing grid counts as spacing decreases. The
manifest must also carry nonempty `caseManifestPaths` and
`caseManifestFilesVerified = true` so each grid case has a verified
fixed-geometry manifest.
`isQuantitativeBenchmark` alone is not sufficient to pass the grid-convergence
requirement.
`precip_AuditYoonBenchmarkReadiness` provides the final conservative gate for
quantitative claims. Current smoke outputs are expected to fail this audit
because they lack complete quantitative reference/geometry packages,
production grid-convergence evidence, and production-validated Stokes
provenance. Reference evidence must include a complete digitization package
with source figure, screenshot, WebPlotDigitizer project, calibration, raw
export, conversion script, uncertainty record, converted CSV, note, and empty
`missingAssets`. The package must also carry `assetFilesVerified = true` from
the loader, verified calibration/raw-export/uncertainty tables, and
`numReferenceRows > 0`, so a hand-authored provenance shell cannot pass
without verified files, machine-readable digitization tables, and nonempty
converted reference data.
`isQuantitativeBenchmark = true` alone is not sufficient.
Fixed-geometry evidence must also carry
`flowFeedbackAccepted = true`, which requires an observed `Vm >= 0.6`
topology change, a recomputed Stokes flow field, and
`flowBackendProductionValidated = true` derived from the final flow evidence.
The fixed-geometry runner now computes that validation state from the
Stokes/non-proxy solver provenance, nonempty `finalFlowSolver`,
`finalFlowLinearResidualRelative`, `maxAcceptedFlowLinearResidualRelative`, and
the recorded residual threshold instead of accepting a hand-authored option.
It also records and gates `finalFlowMaxDivergenceResidual_s_inv` against
`maxAcceptedFlowDivergenceResidual_s_inv`, and
`finalFlowBoundaryClosureRelativeError` against
`maxAcceptedFlowBoundaryClosureRelativeError`.
The readiness gate checks those fields directly, including topology-change
count and recomputation count. It also requires finite post-feedback flow
metrics: `finalRelativePermeability` between 0 and 1,
`finalPressureDropRelative >= 1`, and `finalFlowRateRelative >= 0`. In
addition, `production_stokes` now requires the independent
`yoon_production_flow_validation_manifest.json` to accept the required HyPHM
Stokes validation cases before final benchmark claims. Geometry
evidence must also include nonempty package directory/name/note, substrate mask
file, and region mask file fields from the calibrated geometry package. It must
also include `geometryPackageAssetFilesVerified = true`,
`geometryPackageNumSubstrateCells > 0`, `geometryPackageNumFirstPoreCells > 0`,
and `geometryPackageNumFirstThreePoresCells >=
geometryPackageNumFirstPoreCells`; a standalone
`geometryPackageIsQuantitative = true` flag is not sufficient. Final
reaction-mass evidence must also include finite Ca, C, and alkalinity relative
closure errors no greater than `1e-10`; `finalReactionMassAccepted = true`
alone is not sufficient. Center-band morphology evidence must include positive
active precipitate area and finite centroid-distance/tolerance fields, with the
final precipitate centroid no farther from the split inlet than the recorded
tolerance. Rate-constant evidence must include finite nonnegative
`yoonRateK1`, `yoonRateK2`, and `yoonRateK3` values with
`yoonRateUnits = mol_cm-2_s-1`, nonempty `yoonRateSourceDoi`, nonempty
`yoonRateSourceEquation`, `yoonRateMineralPhase = Vaterite`, and
`yoonRateSourceValuesVerified = true`; a clean-looking source string alone is
not enough. Evidence that used smoke/default/pending Vaterite rate constants or
a finite concentration limiter
also fails the audit until literature-locked rates and a no-clipping production
run are available. Chemistry evidence must carry
`chemistrySpeciationAccepted = true`, `hasIphreeqc = true`, and
`isQuantitativeAcceptance = true`; the readiness gate also checks finite
`maxAbsPhDifference` and `maxAbsSiVateriteDifference` values against their
recorded finite acceptance thresholds. It additionally requires PHREEQC inlet
pH endpoints to match the expected pH 6.1 and 10.9 within the recorded inlet-pH
tolerance, and requires both Yoon and PHREEQC Vaterite SI/omega peaks to occur
at an internal mixing fraction rather than at either inlet endpoint.
Passive transport evidence must carry `passiveTransportAccepted = true` from
the no-reaction split-inlet conservative-transport manifest. The readiness gate
also checks the underlying evidence directly: no finite concentration limiter,
finite nonnegative component concentrations, inert Na/Cl mass error within its
recorded threshold, mixing symmetry error within its threshold, and rejected
transport substeps no greater than the accepted limit. The manifest must also
record the conservative component list `Ca_total`, `C_total`, `Na_total`,
`Cl_total`, and `Alkalinity`, and must not include `H_total` or free-H transport
fields.
Case 1/Case 5 comparison evidence must carry `case1Case5Accepted = true`; the
current fixed-geometry comparison manifest is still marked false because it is
a diagnostic scaffold rather than the production 13/18/118 min comparison. Its
criteria fields make that failure auditable instead of opaque. The readiness
gate checks those criteria fields directly: target-time completeness,
reaction-mass acceptance, delayed/dissolution-only Case 5 factor provenance,
Case 5 boost activation/time evidence, final precipitate/area trend checks,
finite nonnegative Case 1/Case 5 final
moles and area values supporting those trends, and
`productionComparisonValidated` must all be present and accepted.
Diffusion-feedback sensitivity evidence must carry
`diffusionFeedbackAccepted = true` from the Case 2/1/3 manifest. The readiness
gate also checks the underlying evidence directly: `diffusionExponent` must be
`[0, 2, 3]`, `areaTrendNonincreasing` must be true, and
`massTrendDecreasing` must be true. The manifest must also carry finite
nonnegative three-case arrays for `totalPrecipitatedArea_cm2`,
`totalPrecipitateMoles`, and `finalMeanEffectiveDiffusivity_cm2_s`, plus
`feedbackCoupled = true` for all three cases, and the numeric area/mass arrays
must support the recorded trends.
The quantitative-manifest gate requires `isQuantitativeBenchmark = true` plus
nonempty fixed-geometry output evidence fields: `outputRoot`, `snapshotDir`,
`areaCsv`, `flowDiagnosticsCsv`, `transportDiagnosticsCsv`,
`reactionMassLedgerCsv`, `reactionDiagnosticsCsv`, and a nonempty `matFiles`
snapshot list. `outputEvidenceFilesVerified = true` must also confirm those
paths were checked after export. It also requires `transportMode =
'finite_volume'` and `reactionSubcycling = true`. A standalone quantitative
flag is not sufficient.
`precip_RunYoonBenchmarkReadinessAudit` writes the same gate result to CSV/JSON
so each benchmark run can carry an explicit machine-readable pass/fail record.

## Known Limitations

- The current geometry is an approximate local cylindrical-post layout, not a
  digitized Zhang/Yoon micromodel geometry.
- `reference_data/zhang_yoon_reference_curves.csv` contains approximate visual
  digitizations. Use higher-precision WebPlotDigitizer/source-table data before
  making quantitative claims.
- First-pore and first-three-pore areas default to local x-window diagnostics
  unless explicit `spec.regionMasks` are supplied. Exact masks from the
  published figures still need to be digitized and versioned.
- Finite `phreeqcTransportMaxFactor` values are numerical guards against
  transport overshoot before PHREEQC calls. They are not Zhang/Yoon parameters,
  and outputs that depend on them are rejected by the readiness gate for
  quantitative benchmark claims.
- The original dissolution solver has unrelated worktree modifications in the
  current repository state; precipitation development should remain isolated in
  this folder.
- The Yoon micro-continuum path still needs literature-locked rate constants,
  production-scale verification of the CFL time-step policy, a
  production-validated Stokes backend for the `flowSolverFcn` hook after
  `Vm >= 0.6` topology changes, exact area masks, and production 10/5/2.5 um
  grid-convergence runs before quantitative benchmark claims.

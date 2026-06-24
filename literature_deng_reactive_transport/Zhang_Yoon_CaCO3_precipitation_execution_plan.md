# Zhang-Yoon CaCO3 Precipitation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Every key implementation task should be followed by two reviews: parameter/source-basis review and reviewer-perspective code-risk review.

**Goal:** Build a CaCO3 precipitation simulation path in RTSPHEM that reuses the existing dissolution PNM framework while keeping precipitation code isolated under `ReactiveTransport/RTM/precipitate`.

**Architecture:** Copy the existing `PNM_beauty3.m` into a precipitation-specific solver and make all precipitation-specific changes in that copy. Couple a signed PHREEQC calcite reaction interface to level-set motion, add split-inlet CaCl2/Na2CO3 transport, run a local Zhang/Yoon micromodel benchmark, and export precipitation-area/time comparisons against Zhang experiment and Yoon Case 1/Case 5.

**Tech Stack:** MATLAB R2025b, RTSPHEM/HyPHM transport and level-set code, IPhreeqcCOM/PHREEQC, existing RTM helper functions, Markdown documentation.

---

## Current Status Summary

| Step | Status | Evidence |
|---|---|---|
| Copy `PNM_beauty3.m` to `precip_PNM_beauty3.m` | Done | `ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m` exists; first line is `function result = precip_PNM_beauty3(config)` |
| Ensure runner calls copied solver | Done | `run_zhang_yoon_caco3_precipitation_benchmark.m` calls `precip_PNM_beauty3(cfg)` |
| Create Zhang/Yoon runner and config | Done | `run_zhang_yoon_caco3_precipitation_benchmark.m`; `precip_ConfigureZhangYoonBenchmark.m` |
| Implement signed PHREEQC single-cell helpers | Done for helper/integration scope | signed parser, scaler, prescribed TST-match refill, interface-rate helper, signed batch runner, and COM directional tests exist |
| Connect signed PHREEQC into copied PNM | Smoke verified | `useSignedCalciteSurface` selects `precip_RunPhreeqcCalciteBatchSigned`; short Zhang/Yoon run generated `phreeqc_results/phreeqc_summary_log.csv` |
| Level-set sign tests | Done for helper, LevelSetSolver2ndOrder smoke, and short PNM smoke scope | helper-level sign unit test, real LevelSetSolver2ndOrder geometry directional smoke, and short copied-PNM signed-PHREEQC smoke are complete; full benchmark-scale signed coupling remains Task 9 |
| Implement split inlet transport | Done in copied PNM | `precipBuildBoundaryFluxFunction` uses `inletA.yRange` / `inletB.yRange`; mesh-alignment guard has direct tests |
| Run local Zhang/Yoon geometry | Diagnostic full-output run completed; stable benchmark pending | clipped diagnostic PHREEQC run `zhang_yoon_caco3_25mM_20260624_095000` reached 7080 s but stability flags prevent quantitative validation |
| Export precipitation area-time curve and 13/18/118 min images | Done for clipped diagnostic run; stable visual audit pending | `benchmark_snapshot_013min.png`, `benchmark_snapshot_018min.png`, `benchmark_snapshot_118min.png`, `precipitation_area_timeseries.csv/png` exist in `..._095000` |
| Compare with Zhang experiment and Yoon Case 1/Case 5 | Done for clipped diagnostic run; stable scientific comparison pending | `zhang_yoon_area_comparison.csv/png` and `benchmark_comparison_report.md` exist in `..._095000`; report marks the run diagnostic because of finite limiter and stability flags |

## User-Requested Step Checklist

| User-requested step | Current status | Notes / evidence |
|---|---|---|
| Copy `ReactiveTransport/RTM/PNM_beauty3.m` to `ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m` | Done | Precipitation-specific copy exists and exposes `precip_PNM_beauty3(config)` |
| Rename main function to `precip_PNM_beauty3` and make runner call the copy | Done | `run_zhang_yoon_caco3_precipitation_benchmark.m` calls `precip_PNM_beauty3(cfg)` |
| Create `run_zhang_yoon_caco3_precipitation_benchmark.m` | Done | Main benchmark entry point exists under `ReactiveTransport/RTM/precipitate` |
| Create `precip_ConfigureZhangYoonBenchmark.m` | Done | Local Zhang/Yoon benchmark config exists, including split inlet chemistry and benchmark output settings |
| Implement `precip_RunPhreeqcCalciteBatchSigned.m` for single-cell signed PHREEQC | Done for helper/integration scope | Signed PHREEQC runner exists and is covered by MATLAB tests; it now calls the precip-local input builder |
| Create `precip_ComputeSignedCalciteInterfaceRatePerArea.m` | Done | Signed delta-to-interface-rate helper exists and is tested |
| Create `precip_ScaleSignedCalciteDeltaToCellInventory.m` | Done | Signed per-cell scaling helper exists and is tested |
| Connect signed PHREEQC inside `precip_PNM_beauty3.m` | Done for smoke/diagnostic path | `mineralEvolutionMode = 'signed_calcite_surface'` / `useSignedCalciteSurface` selects the signed PHREEQC path |
| Add level-set positive/negative velocity tests | Done for helper and real level-set smoke | Helper tests and `precip_RunLevelSetSignedVelocitySmoke.m` verify precipitation/dissolution sign direction |
| Implement split inlet transport, with optional `precip_CreateTransportMultiInlet.m` | Done | `precip_CreateTransportMultiInlet.m` exists; copied PNM uses split CaCl2/Na2CO3 inlet ranges |
| Run local Zhang/Yoon geometry | Diagnostic run done; stable validation pending | Clipped full-output run reached 118 min; no-clipping hard-CFL run still fails from transport/PHREEQC instability |
| Export precipitation area-time curve and 13/18/118 min images | Done for diagnostic run; stable outputs pending | CSV/PNG area series and benchmark snapshots exist for the clipped diagnostic run |
| Compare with Zhang experiment and Yoon Case 1/Case 5 | Diagnostic comparison done; quantitative validation pending | Comparison CSV/PNG/Markdown report exist, but explicitly disclose limiter/stability flags |

Verification already completed:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate','-begin'); addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\couplePhreeqc'); results = runtests('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate\tests\test_precip_signed_helpers.m'); assertSuccess(results); fprintf('tests=%d passed=%d\n', numel(results), nnz([results.Passed]));"
```

Observed result after split-inlet, PHREEQC single-cell, level-set smoke, mesh-alignment, precipitation-area helper, area time-series export, benchmark snapshot export, comparison-scaffold, `entire_domain` comparison support, source/case/region guard, Zhang selected-pore schema support, repository reference CSV semantics, header-only reference guard, legacy area-CSV compatibility, PHREEQC transport upper-bound helper, adaptive time-grid extension helper, PHREEQC species transport extension helper, precip-local PHREEQC input builder, precip-local PHREEQC output helpers, split-inlet stability-diagnostic inlet flux helper, split-inlet concentration-CFL scale helper, and comparison-report numerical diagnostics additions: 69 MATLAB helper/integration tests passed on the local machine.

Not yet completed:

- Physically stable local benchmark PHREEQC run to 13/18/118 min. A clipped diagnostic full-output run reached 118 min, but it used `phreeqcTransportMaxFactor = 2`, `enableConcentrationCflLimit = false`, and had persistent `overshoot_c;advective_cfl_gt_1;mass_balance_drift` flags.
- Quantitative 13/18/118 min benchmark validation.
- Physically meaningful full-run comparison against the digitized Zhang/Yoon reference rows; current comparison outputs are diagnostic only, and the reference rows are approximate visual digitizations that still need WebPlotDigitizer/source-table refinement before quantitative claims.
- Real figure/image comparison against Zhang and Yoon reference cases.
- Original `ReactiveTransport/RTM/PNM_beauty3.m` currently has unrelated/unreviewed worktree modifications; before merging the precipitation module, audit whether those changes belong to a separate task and keep the precipitation path isolated in `ReactiveTransport/RTM/precipitate`.

---

## File Map

### Existing source used as baseline

- `ReactiveTransport/RTM/PNM_beauty3.m`  
  Original dissolution-oriented PNM solver. This file must not be edited for the precipitation module.

- `ReactiveTransport/RTM/couplePhreeqc/BuildCalcitePhreeqcInput.m`  
  Original shared PHREEQC input builder used as the implementation reference. The signed precipitation path now calls the precip-local copy `precip_BuildCalcitePhreeqcInputSigned.m` instead of depending directly on this shared helper.

### Precipitation module files

- `ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m`  
  Copied precipitation-specific main solver. All split-inlet, signed PHREEQC, and Zhang/Yoon output logic belongs here or in `precip_` helpers.

- `ReactiveTransport/RTM/precipitate/run_zhang_yoon_caco3_precipitation_benchmark.m`  
  User-facing benchmark runner. It builds the config and calls `precip_PNM_beauty3(cfg)`.

- `ReactiveTransport/RTM/precipitate/precip_ConfigureZhangYoonBenchmark.m`  
  Benchmark configuration: geometry, flow, split inlet chemistry, PHREEQC settings, comparison times, and output paths.

- `ReactiveTransport/RTM/precipitate/precip_RunPhreeqcCalciteBatchSigned.m`  
  Signed PHREEQC batch runner that preserves positive precipitation and negative dissolution deltas.

- `ReactiveTransport/RTM/precipitate/precip_BuildCalcitePhreeqcInputSigned.m`  
  Precipitation-local PHREEQC input builder for signed calcite runs. It isolates the signed precipitation path from the modified shared `couplePhreeqc/BuildCalcitePhreeqcInput.m` helper while preserving the selected-output `KIN_DELTA("Calcite")` headings needed by the signed parser.

- `ReactiveTransport/RTM/precipitate/precip_WritePhreeqcSpeciesTable.m`  
  Precipitation-local PHREEQC species CSV writer. It preserves the shared output columns and adds signed calcite diagnostics such as `calciteDeltaMoles`, `calcitePrecipitatedMoles`, and `calciteSignedRate_mol_s`.

- `ReactiveTransport/RTM/precipitate/precip_ExportPhreeqcSpeciesPlots.m`  
  Precipitation-local PHREEQC diagnostic plot exporter for concentration fields and signed calcite amount/rate fields.

- `ReactiveTransport/RTM/precipitate/precip_PrepareConcentrationFaceData.m`  
  Precipitation-local P0 triangle mask helper used by PHREEQC diagnostic plots. It keeps the original all-pore display convention and masks mixed/solid triangles.

- `ReactiveTransport/RTM/precipitate/precip_ParsePhreeqcSelectedOutputSigned.m`  
  Parser for PHREEQC selected output that keeps `KIN_DELTA("Calcite")` sign.

- `ReactiveTransport/RTM/precipitate/precip_ScaleSignedCalciteDeltaToCellInventory.m`  
  Scales PHREEQC reference-water reaction amounts to per-cell inventory and caps impossible dissolution.

- `ReactiveTransport/RTM/precipitate/precip_ApplyPrescribedCalciteDissolutionSigned.m`  
  Restores prescribed TST-match dissolution amounts when PHREEQC input uses prescribed `REACTION` rather than kinetic `KIN_DELTA`.

- `ReactiveTransport/RTM/precipitate/precip_ComputeSignedCalciteInterfaceRatePerArea.m`  
  Converts signed mineral amount into level-set interface velocity.

- `ReactiveTransport/RTM/precipitate/precip_ComputePhreeqcTransportUpperBounds.m`  
  Computes optional pre-PHREEQC transport concentration caps from initial fields, legacy scalar inlets, and Zhang/Yoon split-inlet chemistry. A finite cap is a numerical guard against transport overshoot, not a literature parameter.

- `ReactiveTransport/RTM/precipitate/precip_ResolveAdaptiveTimeGridExtension.m`  
  Strictly parses the optional `allowAdaptiveTimeGridExtension` runtime control. This is a numerical time-stepper control, not a Zhang/Yoon physical parameter.

- `ReactiveTransport/RTM/precipitate/precip_ExtendTransportProblemsToStepper.m`  
  Extends all active transport variables after appending an adaptive internal time step, including the H/Ca/C/Na/Cl transport problems used by signed PHREEQC runs.

- `ReactiveTransport/RTM/precipitate/precip_ComputeBoundaryInletMassFlux.m`  
  Computes the diagnostic inlet solute flux used by stability mass-balance checks. Uniform mode preserves the legacy scalar formula; split-left inlet mode integrates `inletA` and `inletB` segment concentrations over their configured `yRange` lengths.

- `ReactiveTransport/RTM/precipitate/precip_ComputeBoundaryInletConcentrationScale.m`  
  Computes the diagnostic inlet concentration scale used by `negative_c` and `overshoot_c` stability flags. Split-left inlet mode uses the maximum configured inlet concentration for the requested species.

- `ReactiveTransport/RTM/precipitate/precip_CreateTransportMultiInlet.m`  
  Split-left inlet boundary helper for CaCl2 and Na2CO3 streams, with validation for species, concentration values, y-ranges, and overlap.

- `ReactiveTransport/RTM/precipitate/precip_ValidateSplitInletMeshAlignment.m`  
  Verifies that split-inlet y-range endpoints coincide with left-boundary mesh nodes so one boundary segment is not assigned two inlet chemistries.

- `ReactiveTransport/RTM/precipitate/precip_RunLevelSetSignedVelocitySmoke.m`  
  Real `LevelSetSolver2ndOrder` smoke helper that verifies negative normal velocity grows solid and positive velocity shrinks solid in the precipitation sign convention.

- `ReactiveTransport/RTM/precipitate/precip_ComputePrecipitationAreaMetrics.m`  
  Computes total, first-pore, and first-three-pore solid/net-solid area metrics. First-pore windows are approximate barycenter x-window diagnostics, not digitized Zhang/Yoon masks.

- `ReactiveTransport/RTM/precipitate/precip_BenchmarkSnapshotFilename.m`  
  Centralizes stable benchmark snapshot filename rules, including `013min`, `018min`, `118min`, integer-second, and decimal-second artificial smoke names.

- `ReactiveTransport/RTM/precipitate/precip_ExportBenchmarkSnapshots.m`  
  Exports benchmark snapshot PNGs with current interface, initial interface, lower/upper split-inlet labels, simulation time, and stable filenames.

- `ReactiveTransport/RTM/precipitate/precip_CompareZhangYoonBenchmark.m`  
  Compares simulated precipitation-area time series with Zhang/Yoon reference curves, writes comparison CSV/PNG/Markdown report, refuses to run on a header-only reference CSV, and reports finite transport limiters plus stability diagnostic flags to prevent overclaiming.

- `ReactiveTransport/RTM/precipitate/reference_data/zhang_yoon_reference_curves.csv`  
  Source-bounded reference-curve schema with 18 approximate visual digitization rows for 13/18/118 min. Zhang rows store normalized pixel-area values; Yoon rows store area converted from square micrometers to square centimeters.

- `ReactiveTransport/RTM/precipitate/reference_data/digitization_notes.md`  
  Documents the source figures, approximate axis calibration, unit conversion, and accuracy limits for the repository reference rows.

- `ReactiveTransport/RTM/precipitate/reference_data/README.md`  
  Documents that placeholder or synthetic values must not be added to the repository reference CSV.

- `ReactiveTransport/RTM/precipitate/tests/test_precip_signed_helpers.m`  
  Lightweight MATLAB tests for config compatibility, signed parser/scaler, prescribed reaction refill, interface-rate sign, area metrics/time-series export, benchmark snapshot export, and comparison-script input/output guards.

- `ReactiveTransport/RTM/precipitate/README.md`  
  Module handoff documentation for entry points, signed calcite interpretation, outputs, current diagnostic-run status, no-clipping hard-CFL failure evidence, and known limitations.

### Files/data still to complete

- Higher-precision reference digitization  
  The repository CSV now has approximate visual digitization rows for the required comparison times. Before making quantitative scientific claims, refine these rows with WebPlotDigitizer or publisher/source table data and update `digitization_notes.md`.

- Optional: `ReactiveTransport/RTM/precipitate/README.md`  
  Module-level user documentation for running tests, smoke runs, and full benchmark runs.

---

## Task 1: Copy PNM Solver and Lock Module Boundary

**Files:**

- Source: `ReactiveTransport/RTM/PNM_beauty3.m`
- Create/modify: `ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m`
- Create/modify: `ReactiveTransport/RTM/precipitate/run_zhang_yoon_caco3_precipitation_benchmark.m`

- [x] **Step 1: Copy original solver**

Create:

```text
ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m
```

from:

```text
ReactiveTransport/RTM/PNM_beauty3.m
```

Do not edit `ReactiveTransport/RTM/PNM_beauty3.m`.

- [x] **Step 2: Rename copied solver function**

The first line must be:

```matlab
function result = precip_PNM_beauty3(config)
```

- [x] **Step 3: Create benchmark runner**

The runner must end with:

```matlab
cfg = precip_ConfigureZhangYoonBenchmark(overrides);
result = precip_PNM_beauty3(cfg);
```

It must not call:

```matlab
PNM_beauty3(cfg)
```

- [x] **Step 4: Verify function resolution**

Run:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate'); which precip_PNM_beauty3; which run_zhang_yoon_caco3_precipitation_benchmark"
```

Expected:

```text
...ReactiveTransport\RTM\precipitate\precip_PNM_beauty3.m
...ReactiveTransport\RTM\precipitate\run_zhang_yoon_caco3_precipitation_benchmark.m
```

- [x] **Step 5: Review**

Parameter/source review: no literature parameter is introduced in this task.  
Reviewer-risk review: confirm no original solver code was changed for precipitation.

Current status: completed.

---

## Task 2: Build Zhang/Yoon Benchmark Config

**Files:**

- Create/modify: `ReactiveTransport/RTM/precipitate/precip_ConfigureZhangYoonBenchmark.m`
- Test: `ReactiveTransport/RTM/precipitate/tests/test_precip_signed_helpers.m`

- [x] **Step 1: Define core benchmark fields**

The config must include:

```matlab
cfg.benchmarkName = 'zhang_yoon_caco3_25mM';
cfg.reactionModel = 'phreeqc';
cfg.phreeqcRunGroup = 'phreeqc_database_calcite';
cfg.phreeqcRateLaw = 'phreeqc_database_calcite';
cfg.mineralEvolutionMode = 'signed_calcite_surface';
cfg.layoutType = 'zhang2010_micromodel_local';
```

- [x] **Step 2: Define literature-based geometry and transport parameters**

Use:

```matlab
cfg.lengthXAxis = 0.30;
cfg.lengthYAxis = 0.30;
cfg.thickness = 0.002;
cfg.thicknessCm = cfg.thickness;
cfg.postDiameter = 0.03;
cfg.poreBody = 0.018;
cfg.poreThroat = 0.004;
cfg.targetPorosity = 0.39;
cfg.inletVelocity = 0.020833333;
cfg.diffusionCoefficient = 0.9e-5;
```

Parameter basis:

- `inletVelocity = 1.25 cm/min = 0.020833333 cm/s`, Zhang/Yoon benchmark.
- `diffusionCoefficient = 0.9e-5 cm2/s`, Yoon simulation value.
- `thickness = 0.002 cm`, approximation of 17-21 um micromodel depth.
- `postDiameter = 0.03 cm`, 300 um posts.
- `poreBody = 0.018 cm`, 180 um pore body.
- `poreThroat = 0.004 cm`, 40 um throat.
- `lengthXAxis = lengthYAxis = 0.30 cm` is a local crop, not the full Zhang 2 cm x 1 cm domain.

- [x] **Step 3: Define split-inlet chemistry**

Use:

```matlab
cfg.flowBoundaryMode = 'split_left_inlet';
cfg.splitInletY = 0.5 * cfg.lengthYAxis;

cfg.inletA.name = 'CaCl2';
cfg.inletA.yRange = [0, cfg.splitInletY];
cfg.inletA.H_total = 10^(-6.1) * 1e-3;
cfg.inletA.Ca_total = 25e-6;
cfg.inletA.C_total = 0;
cfg.inletA.Na_total = 0;
cfg.inletA.Cl_total = 50e-6;

cfg.inletB.name = 'Na2CO3';
cfg.inletB.yRange = [cfg.splitInletY, cfg.lengthYAxis];
cfg.inletB.H_total = 10^(-10.9) * 1e-3;
cfg.inletB.Ca_total = 0;
cfg.inletB.C_total = 25e-6;
cfg.inletB.Na_total = 50e-6;
cfg.inletB.Cl_total = 0;
```

Unit conversion:

```text
25 mM = 25e-6 mol/cm3
50 mM = 50e-6 mol/cm3
```

- [x] **Step 4: Preserve legacy scalar inlet fields for compatibility**

Set premixed fallback fields so old scalar code paths do not silently use zero:

```matlab
cfg.inletCalciumConcentration = 0.5 * cfg.inletA.Ca_total;
cfg.inletCarbonConcentration = 0.5 * cfg.inletB.C_total;
cfg.inletSodiumConcentration = 0.5 * cfg.inletB.Na_total;
cfg.inletChlorideConcentration = 0.5 * cfg.inletA.Cl_total;
```

- [x] **Step 5: Add config compatibility test**

Add test:

```matlab
function testConfigureBenchmarkKeepsPnmCompatibilityFields(testCase)
cfg = precip_ConfigureZhangYoonBenchmark(struct('lengthYAxis', 0.4, 'thickness', 0.003));

verifyEqual(testCase, cfg.thicknessCm, 0.003, 'AbsTol', 0);
verifyEqual(testCase, cfg.splitInletY, 0.2, 'AbsTol', 0);
verifyEqual(testCase, cfg.targetLengthYAxis, cfg.lengthYAxis, 'AbsTol', 0);
verifyEqual(testCase, cfg.targetAspectRatio, cfg.lengthXAxis / cfg.lengthYAxis, 'AbsTol', 0);
verifyEqual(testCase, cfg.inletCalciumConcentration, 12.5e-6, 'AbsTol', 1e-15);
verifyEqual(testCase, cfg.inletCarbonConcentration, 12.5e-6, 'AbsTol', 1e-15);
verifyEqual(testCase, cfg.mineralEvolutionMode, 'signed_calcite_surface');
end
```

- [x] **Step 6: Review**

Parameter/source review: confirm all values are tied to Zhang/Yoon or marked as local-crop approximation.  
Reviewer-risk review: confirm config fields match copied PNM field names.

Current status: completed.

---

## Task 3: Implement Single-Cell Signed PHREEQC Helpers

**Files:**

- Create/modify: `ReactiveTransport/RTM/precipitate/precip_RunPhreeqcCalciteBatchSigned.m`
- Create/modify: `ReactiveTransport/RTM/precipitate/precip_ParsePhreeqcSelectedOutputSigned.m`
- Create/modify: `ReactiveTransport/RTM/precipitate/precip_ScaleSignedCalciteDeltaToCellInventory.m`
- Create/modify: `ReactiveTransport/RTM/precipitate/precip_ApplyPrescribedCalciteDissolutionSigned.m`
- Test: `ReactiveTransport/RTM/precipitate/tests/test_precip_signed_helpers.m`

- [x] **Step 1: Parse signed PHREEQC output**

`precip_ParsePhreeqcSelectedOutputSigned` must preserve:

```matlab
result.calciteDeltaMoles = KIN_DELTA_Calcite;
result.calcitePrecipitatedMoles = max(result.calciteDeltaMoles, 0);
result.calciteDissolvedMoles = max(-result.calciteDeltaMoles, 0);
result.calciteSignedRate_mol_s = result.calciteDeltaMoles ./ timeStepSize;
```

- [x] **Step 2: Scale signed delta to cell inventory**

`precip_ScaleSignedCalciteDeltaToCellInventory` must:

- Scale PHREEQC reference-water delta by cell water volume.
- Prevent dissolution from exceeding available `calcite_moles`.
- Keep precipitation uncapped by initial aqueous inventory by default, because PHREEQC already returns post-reaction water chemistry.
- If optional precipitation capping is enabled, blend water chemistry fields from real pre-state fields rather than zero.

- [x] **Step 3: Restore prescribed TST-match dissolution**

`precip_ApplyPrescribedCalciteDissolutionSigned` must set:

```matlab
result.calciteDeltaMoles = -dissolvedMoles;
result.calcitePrecipitatedMoles = zeros(numCells, 1);
result.calciteDissolvedMoles = dissolvedMoles;
result.calciteSignedRate_mol_s = -dissolvedMoles ./ max(timeStepSize, eps);
result.calciteRate_mol_s = dissolvedMoles ./ max(timeStepSize, eps);
```

- [x] **Step 4: Run single-cell helper tests without PHREEQC COM**

Run:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate'); addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\couplePhreeqc'); results = runtests('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate\tests\test_precip_signed_helpers.m'); assertSuccess(results);"
```

Expected:

```text
test_precip_signed_helpers
.................
```

- [x] **Step 5: Run single-cell PHREEQC COM test**

Original conceptual half-concentration smoke state:

```matlab
state = struct();
state.h_mol_cm3 = 10^(-8.3) * 1e-3;
state.ca_mol_cm3 = 12.5e-6;
state.c_mol_cm3 = 12.5e-6;
state.na_mol_cm3 = 25e-6;
state.cl_mol_cm3 = 25e-6;
state.interface_area_cm2 = 1e-3;
state.water_volume_cm3 = 1e-3;
state.calcite_moles = 1e-6;

options = struct();
options.databasePath = 'C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat';
options.workDir = tempname;
options.timeStepIndex = 1;
options.timeStepSize = 1;
options.rateLaw = 'database_calcite';

result = precip_RunPhreeqcCalciteBatchSigned(state, options);
disp(result.calciteDeltaMoles);
disp(result.calciteSI);
```

Expected:

- Function completes without COM/database error.
- `result.calciteDeltaMoles` is finite.
- `result.pH`, `result.ca_total_mol_cm3`, `result.c_total_mol_cm3`, and `result.calciteSI` are finite for the active cell.
- Local observed smoke result with USGS `phreeqc.dat`: `calciteDeltaMoles = 0`, `calciteSI = 0.8281`.

Automated test:

```matlab
testSignedPhreeqcSingleCellLocalDatabaseSmoke
```

Scope:

- The conceptual premixed single-cell state above is a representative half-concentration smoke state for the Zhang/Yoon 25 mM two-inlet chemistry.
- The automated COM test uses artificial enlarged PHREEQC coupling scales to make signed `KIN_DELTA` observable in a fast unit/integration test:

```matlab
state.interface_area_cm2 = 1e4;
state.water_volume_cm3 = 1000;
state.calcite_moles = 10;
options.timeStepSize = 1000;
options.maxSpecificSurfaceArea = 1e12;
options.kineticsReservoirMoles = 10;
```

- The directional stress cases are not Zhang/Yoon literature parameter values; they are non-literature PHREEQC sign tests for the coupling code.
- The repository-local `phreeqc-m.dat` is absent in the current worktree; local verification used USGS IPhreeqcCOM `phreeqc.dat`.

Directional signed-delta coverage:

```text
supersaturated single cell -> calciteDeltaMoles > 0
undersaturated single cell -> calciteDeltaMoles < 0
```

The test also verifies `USER_PUNCH` / `KIN_DELTA_Calcite RATE_Calcite` input headings, unique `tempname` work directories, `calciteSignedRate_mol_s = calciteDeltaMoles / timeStepSize`, and signed split into `calcitePrecipitatedMoles` / `calciteDissolvedMoles`.

- [x] **Step 6: Review**

Parameter/source review: verify signed `KIN_DELTA("Calcite")` interpretation follows PHREEQC output and existing `BuildCalcitePhreeqcInput` headings.  
Reviewer-risk review: verify no dissolution-only `max(-delta, 0)` path is used to drive precipitation.

Current status: helper tests and real PHREEQC COM single-cell run completed on the local machine. The repository-local `phreeqc-m.dat` was not present, so the verified local environment used USGS IPhreeqcCOM `phreeqc.dat`.

---

## Task 4: Convert Signed Delta to Level-Set Interface Rate

**Files:**

- Create/modify: `ReactiveTransport/RTM/precipitate/precip_ComputeSignedCalciteInterfaceRatePerArea.m`
- Modify: `ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m`
- Test: `ReactiveTransport/RTM/precipitate/tests/test_precip_signed_helpers.m`

- [x] **Step 1: Implement signed interface-rate conversion**

Use:

```matlab
ratePerArea = -deltaMoles ./ max(timeStepSize, eps) ./ interfaceAreaCm2;
```

Sign convention:

```text
calciteDeltaMoles > 0 -> precipitation -> ratePerArea < 0 -> solid grows
calciteDeltaMoles < 0 -> dissolution -> ratePerArea > 0 -> solid shrinks
```

- [x] **Step 2: Add sign unit test**

Test:

```matlab
function testSignedInterfaceRateMakesPrecipitationNegative(testCase)
result = struct();
result.calciteDeltaMoles = [2e-9; -3e-9; 1e-9];
rate = precip_ComputeSignedCalciteInterfaceRatePerArea(result, [1e-3; 2e-3; 0], 10);

verifyEqual(testCase, rate, [-2e-7; 1.5e-7; 0], 'AbsTol', 1e-18);
end
```

- [x] **Step 3: Connect signed rate in copied PNM**

In `precip_PNM_beauty3.m`, signed mode must use:

```matlab
if useSignedCalciteSurface
    phreeqcCalciteRateData = precip_ComputeSignedCalciteInterfaceRatePerArea( ...
        phreeqcSpeciesData, phreeqcGeometryState.interface_area_cm2, macroscaleTimeStepSize);
else
    phreeqcCalciteRateData = ComputePhreeqcInterfaceRatePerArea( ...
        phreeqcSpeciesData, phreeqcGeometryState.interface_area_cm2, ...
        macroscaleTimeStepSize, phreeqcOptions, phreeqcState.h_mol_cm3);
end
```

- [x] **Step 4: Preserve negative rates during level-set vertex averaging**

Signed mode must not use positive-only averaging:

```matlab
if strcmpi(phreeqcRateLaw, 'tst_match') && ~useSignedCalciteSurface
    phreeqcRateV = averagePositiveTriangleDataToVertices(...);
else
    phreeqcRateV = averageTriangleDataToVertices(...);
end
```

- [x] **Step 5: Add level-set geometry smoke test**

Create a minimal test helper or script that:

1. Builds a small level-set circle.
2. Applies a uniform negative rate.
3. Advances one step.
4. Verifies solid area increases.
5. Applies a uniform positive rate.
6. Advances one step.
7. Verifies solid area decreases.

Expected:

```text
precipitation step: solid area after > solid area before
dissolution step: solid area after < solid area before
```

Implemented:

- `ReactiveTransport/RTM/precipitate/precip_RunLevelSetSignedVelocitySmoke.m`
- `testSignedLevelSetSmokeMakesNegativeSpeedGrowSolid`
- `testSignedLevelSetSmokeRejectsInvalidOptions`

The smoke test constructs a circular signed-distance level set on a `FoldedCartesianGrid`, advances it with the real `levelSetEquationTimeStep`, and checks:

```text
negative normal speed -> solid area increases
positive normal speed -> solid area decreases
```

Smoke-test parameters:

```matlab
domainLengthCm = 1.0;
initialRadiusCm = 0.20;
speedMagnitudeCmS = 0.02;
timeStepSize = 1.0;
```

These are artificial smoke-test parameters, not Zhang/Yoon literature values and not benchmark-scale rates. Solid area is estimated by node counting for this lightweight test; production benchmark area metrics should use triangle/subcell geometry. This closes the level-set directional smoke requirement, but it does not replace the later short PNM run that must cover the full pipeline:

```text
phreeqcCalciteRateData -> averageTriangleDataToVertices -> molarVolume * phreeqcRateV -> CFL substeps -> geometry/porosity update
```

- [x] **Step 6: Review**

Parameter/source review: no new literature parameter is introduced.  
Reviewer-risk review: confirm signed precipitation cannot be filtered out by a dissolution-only positive-rate path.

Current status: helper-level sign test and real LevelSetSolver2ndOrder geometry directional smoke test completed. A full PNM signed-coupling smoke remains part of Task 6.

---

## Task 5: Implement Split-Inlet Transport

**Files:**

- Modify: `ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m`
- Optional create: `ReactiveTransport/RTM/precipitate/precip_CreateTransportMultiInlet.m`
- Test: `ReactiveTransport/RTM/precipitate/tests/test_precip_signed_helpers.m` or a new `test_precip_split_inlet.m`

- [x] **Step 1: Add split-inlet boundary mode**

In copied PNM:

```matlab
flowBoundaryMode = lower(strrep(strtrim(char(cfgget(config, ...
    'flowBoundaryMode', 'uniform_left_inlet'))), '-', '_'));
```

- [x] **Step 2: Build flux functions for each transported species**

Use local helper:

```matlab
inletFluxHydrogen = precipBuildBoundaryFluxFunction( ...
    flowBoundaryMode, config, 'H_total', initialHydrogenConcentration, inletVelocity, EPS);
inletFluxCalcium = precipBuildBoundaryFluxFunction( ...
    flowBoundaryMode, config, 'Ca_total', inletCalciumConcentration, inletVelocity, EPS);
inletFluxCarbon = precipBuildBoundaryFluxFunction( ...
    flowBoundaryMode, config, 'C_total', inletCarbonConcentration, inletVelocity, EPS);
inletFluxSodium = precipBuildBoundaryFluxFunction( ...
    flowBoundaryMode, config, 'Na_total', inletSodiumConcentration, inletVelocity, EPS);
inletFluxChloride = precipBuildBoundaryFluxFunction( ...
    flowBoundaryMode, config, 'Cl_total', inletChlorideConcentration, inletVelocity, EPS);
```

- [x] **Step 3: Respect `inletA.yRange` and `inletB.yRange`**

The split flux helper must use:

```matlab
aRange = cfgget(inletA, 'yRange', [0, splitInletY]);
bRange = cfgget(inletB, 'yRange', [splitInletY, cfgget(config, 'lengthYAxis', Inf)]);
```

Boundary split rule:

```matlab
lowerInlet = yCoord >= min(aRange) & yCoord < max(aRange);
upperInlet = yCoord >= min(bRange) & yCoord <= max(bRange);
```

The shared boundary point belongs to `inletB` to avoid double-counting.

- [x] **Step 4: Connect flux functions to transport objects**

`createPhreeqcTransport` must accept an optional `inletFluxFunction`:

```matlab
function transportProblem = createPhreeqcTransport(..., inletFluxFunction)
if nargin < 11 || isempty(inletFluxFunction)
    inletFluxFunction = precip_CreateTransportMultiInlet('uniform_left_inlet', ...
        struct(), '', inletConcentration, inletVelocity, EPS);
end
transportProblem.gF.setdata(inletFluxFunction);
```

- [x] **Step 5: Add split-inlet flux unit test**

Add a small test that evaluates flux at three points:

```matlab
x = [0, 0, 0.01; 0.01, 0.20, 0.20];
```

Expected for Ca:

```text
lower-left point: negative Ca flux
upper-left point: zero Ca flux
non-left point: zero flux
```

Expected for C:

```text
lower-left point: zero C flux
upper-left point: negative C flux
non-left point: zero flux
```

If the local helper remains private inside `precip_PNM_beauty3.m`, extract it to:

```text
ReactiveTransport/RTM/precipitate/precip_CreateTransportMultiInlet.m
```

so this test can call it directly.

Implemented:

- `ReactiveTransport/RTM/precipitate/precip_CreateTransportMultiInlet.m`
- `testSplitInletFluxSeparatesCalciumAndCarbonate`
- `testSplitInletFluxCoversHydrogenSodiumAndChloride`
- `testUniformInletFluxKeepsLegacyScalarBehaviorAndNx2Coordinates`
- `testSplitInletRejectsMissingConfigOrSpecies`
- `testSplitInletRejectsInvalidConcentrations`
- `testSplitInletRejectsInvalidRanges`

Additional guards:

- `split_left_inlet` now errors if `inletA`/`inletB` or species fields are missing.
- Split-inlet concentrations must be finite numeric scalars.
- Split y-ranges must be finite two-element vectors with positive width, and `inletA.yRange(2)` must meet `inletB.yRange(1)` without overlap.
- `split_left_inlet` currently requires `flowDirection = left_to_right`.
- Split y-range endpoints must align with left-boundary mesh nodes before benchmark runs.
- CaCl2 lower inlet and Na2CO3 upper inlet are the RTSPHEM local coordinate convention for the Zhang/Yoon two-stream boundary, not an asserted absolute image orientation from the papers.

- [x] **Step 6: Review**

Parameter/source review: verify split inlet reflects Zhang/Yoon CaCl2/Na2CO3 boundary.  
Reviewer-risk review: verify old scalar fallback remains available for smoke tests.

Current status: implementation and direct split-flux tests completed. Parameter/source subagent review passed. Reviewer-risk subagent review was used to add explicit config/species/range errors, flow-direction guard, and split mesh-alignment guard.

Additional mesh-alignment verification added:

- `precip_ValidateSplitInletMeshAlignment.m` is now a reusable helper.
- `precip_PNM_beauty3.m` delegates its split-left boundary mesh check to that helper.
- `testSplitInletMeshAlignmentRejectsMissingLeftBoundaryEndpoint` verifies that an inlet split endpoint missing from the left boundary mesh raises `RTSPHEM:Precipitate:SplitInletNotMeshAligned`.
- `testSplitInletMeshAlignmentAcceptsRequiredEndpoints` verifies that aligned endpoints are accepted without warning.
- `testSplitInletMeshAlignmentAcceptsCoordVGridStruct` and `testSplitInletMeshAlignmentRejectsCoordVGridStructMissingEndpoint` cover the production-style `coordV` grid path.
- `testSplitInletMeshAlignmentAcceptsTwoByNCoordinates` covers the HyPHM-style `2 x N` coordinate orientation.

---

## Task 6: Run Local Zhang/Yoon Geometry

**Files:**

- Modify: `ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m`
- Modify: `ReactiveTransport/RTM/precipitate/precip_ConfigureZhangYoonBenchmark.m`
- Output: `ReactiveTransport/RTM/precipitate/outputs/rtm_runs/<runName>/`

- [x] **Step 1: Add `zhang2010_micromodel_local` and `zhang2010_micromodel_full` layout cases**

The copied PNM must accept:

```matlab
case {'zhang2010_micromodel_local', 'zhang2010_micromodel_full'}
```

Both cases use Zhang/Yoon dimensions from config.

- [x] **Step 2: Generate approximate cylindrical-post geometry**

Use:

```matlab
postDiameter = cfgget(config, 'postDiameter', 0.03);
poreBody = cfgget(config, 'poreBody', 0.018);
poreThroat = cfgget(config, 'poreThroat', 0.004);
circleRadius = 0.5 * postDiameter;
pitchX = postDiameter + poreBody;
pitchY = postDiameter + poreThroat;
```

Mark this as:

```text
approximate cylindrical-post geometry
```

Strict reproduction should use digitized DXF/TIF geometry.

- [x] **Step 3: Run a short geometry smoke test**

Run:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate'); result = run_zhang_yoon_caco3_precipitation_benchmark(struct('endTime', 2, 'maxTotalTimeSteps', 1, 'saveMainPlot', false, 'saveIndividualPlots', false, 'saveRealtimePlot', false, 'exportDXF', false, 'saveInterfaceMask', false, 'writeExcel', false)); disp(result.resultsDir);"
```

Expected:

- Geometry initializes.
- Mesh diagnostics can be written or intentionally disabled.
- No unknown `layoutType` error.
- No call to original `PNM_beauty3`.

Observed local evidence:

```text
Run directory:
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260623_224132

Key files:
run_metadata.json
global_evolution_log.csv
stability_diagnostics_log.csv
mesh_diagnostics/mesh_statistics.csv
phreeqc_results/phreeqc_summary_log.csv
phreeqc_results/phreeqc_species_0001.csv ... phreeqc_species_0006.csv
interface_images/timestep_0001.png ... timestep_0006.png
```

The smoke run used the local Zhang/Yoon geometry branch, initialized the mesh, generated geometry/interface outputs, and wrote metadata under the precipitation module output tree.

- [x] **Step 4: Run a short signed-PHREEQC local smoke test**

Run with a very short time:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate'); result = run_zhang_yoon_caco3_precipitation_benchmark(struct('endTime', 5, 'maxTotalTimeSteps', 2, 'saveMainPlot', false, 'saveIndividualPlots', false, 'saveRealtimePlot', false, 'exportDXF', false, 'saveInterfaceMask', false, 'writeExcel', false)); disp(result.resultsDir);"
```

Expected:

- PHREEQC COM starts.
- `phreeqc_summary_log.csv` exists.
- `calcite_net_delta_moles_mean`, `calcite_precipitated_moles_mean`, and `calcite_dissolved_moles_mean` columns exist.

Observed local evidence:

```text
phreeqc_results/phreeqc_summary_log.csv columns include:
timestep
time_s
pH_mean
calcite_rate_mean_mol_cm2_s
calcite_signed_rate_mol_s_mean
calcite_net_delta_moles_mean
calcite_precipitated_moles_mean
calcite_dissolved_moles_mean
calcite_si_mean
charge_balance_mean
```

The short run produced 6 PHREEQC summary rows ending at `time_s = 1.05` and did not crash. This verifies the full file-level path:

```text
runner -> copied PNM -> split inlet -> signed PHREEQC -> signed rate summary -> level-set pipeline entry -> outputs
```

Important limitations from this smoke run:

- It is not a physically validated Zhang/Yoon comparison run.
- `maxTotalTimeSteps = 1` was recorded in metadata, but the current adaptive/time-grid logic produced 6 output records; this is acceptable for smoke evidence but should be clarified before using the option as a strict one-step limit.
- PNG-NMR post-processing was enabled by inherited config and produced non-fatal mesh-size warnings; future RTM-only smoke runs should pass `enablePNGSimulation=false`, `enableNMRSimulation=false`, and `enableNMRSurrogate=false`.
- Stability diagnostics flagged concentration overshoot, advective CFL greater than 1, and mass-balance drift; a calibrated benchmark run needs smaller time steps or stricter CFL control before curve comparison.
- The mean signed reaction in this short run is very small and all reported mean `calcite_precipitated_moles_mean` values are positive while `calcite_dissolved_moles_mean` remains zero; this demonstrates a precipitation-signed path but not a mature clogging/front evolution.

Second no-NMR RTM-only smoke evidence:

```text
Run directory:
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260623_225255

Command intent:
endTime = 0.02 s
timeStepperType = linear
initialMacroscaleTimeStepSize = 0.01 s
enablePNGSimulation = false
enableNMRSimulation = false
enableNMRSurrogate = false
nmr_method = none
```

Observed evidence:

- `run_metadata.json` records `nmr.enabled = false` and `nmr.png_enabled = false`.
- No `png_nmr_results` directory was created.
- `phreeqc_results/phreeqc_summary_log.csv` was written with signed calcite columns.
- `stability_diagnostics_log.csv` records advective CFL below 1 (`0.421` then `0.211`), so the earlier `advective_cfl_gt_1` warning was removed for this smoke.
- The same stability log still flags `overshoot_c;mass_balance_drift`, so this is not yet a physically acceptable benchmark run.

- [x] **Step 5: Review**

Parameter/source review: confirm local crop and geometry approximation are documented.  
Reviewer-risk review: inspect mesh quality, initial porosity, and whether precipitated region appears near the mixing interface.

Current status: code path exists, one PNG-enabled short signed-PHREEQC smoke run completed, and one no-NMR RTM-only smoke run completed. Parameter/source subagent review passed for the mesh-alignment safeguard and smoke framing. Reviewer-risk subagent initially failed the narrow mesh-alignment change because `coordV` and `2 x N` coordinate paths were not covered; the helper/test fixes were applied, local tests passed, and reviewer-risk re-review passed for the narrow mesh-alignment/no-NMR smoke scope. Remaining transport overshoot and mass-balance concerns are not blockers for Task 6 smoke wiring, but they are open Task 9 full-benchmark risks.

---

## Task 7: Export Precipitation Area-Time Curve and 13/18/118 min Images

**Files:**

- Modify or create: `ReactiveTransport/RTM/precipitate/precip_ComputePrecipitationAreaMetrics.m`
- Modify or create: `ReactiveTransport/RTM/precipitate/precip_ExportBenchmarkSnapshots.m`
- Modify: `ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m`
- Output: `benchmark_comparison_times_log.csv`, `precipitation_area_timeseries.csv`, `precipitation_area_timeseries.png`, `benchmark_snapshots/*.png`

- [x] **Step 1: Consume benchmark comparison times**

The solver reads:

```matlab
benchmarkComparisonTimes_s = cfgget(config, 'benchmarkComparisonTimes_s', []);
```

Config sets:

```matlab
cfg.benchmarkComparisonTimes_s = [13, 18, 118] * 60;
```

- [x] **Step 2: Write benchmark comparison CSV**

Current output:

```text
benchmark_comparison_times_log.csv
```

Columns include:

```text
target_time_s
target_time_min
captured_timestep
captured_time_s
porosity
permeability_mD
k_k0
rate_mol_cm2_s
surface_area_cm2
grain_volume_cm3
first_pore_net_solid_area_cm2
first_three_pores_net_solid_area_cm2
first_pore_solid_area_cm2
first_three_pores_solid_area_cm2
```

- [x] **Step 3: Create reusable precipitation area helper**

Created:

```text
ReactiveTransport/RTM/precipitate/precip_ComputePrecipitationAreaMetrics.m
```

Function signature:

```matlab
function metrics = precip_ComputePrecipitationAreaMetrics(grid, levelSetNow, levelSetInitial, config)
```

Return:

```matlab
metrics = struct();
metrics.totalSolidArea_cm2 = ...;
metrics.totalNetSolidArea_cm2 = ...;
metrics.firstPoreSolidArea_cm2 = ...;
metrics.firstPoreNetSolidArea_cm2 = ...;
metrics.firstThreePoresSolidArea_cm2 = ...;
metrics.firstThreePoresNetSolidArea_cm2 = ...;
```

Additional returned fields document the approximate window definition:

```matlab
metrics.firstPoreXMaxCm = ...;
metrics.firstThreePoresXMaxCm = ...;
metrics.windowDefinition = 'approximate_x_windows';
```

Implemented tests:

- `testPrecipitationAreaMetricsComputesTotalAndWindowNetArea`
- `testPrecipitationAreaMetricsUsesZhangWindowDefaults`
- `testPrecipitationAreaMetricsComputesClippedTriangleArea`
- `testPrecipitationAreaMetricsRejectsMalformedGrid`
- `testPrecipitationAreaMetricsRejectsInvalidLevelSets`

Review result:

- Parameter/source re-review passed because no new literature parameter is introduced; the helper explicitly marks first-pore and first-three-pore windows as approximate barycenter x-window diagnostics for the local cylindrical-post layout.
- Reviewer-risk re-review passed after hardening grid and level-set validation. The helper was accepted as created/tested; subsequent Step 4/5 work wires it into production time-series export.

- [x] **Step 4: Export precipitation area time series**

Create:

```text
precipitation_area_timeseries.csv
```

Required columns:

```text
timestep,time_s,total_net_solid_area_cm2,first_pore_net_solid_area_cm2,first_three_pores_net_solid_area_cm2,total_solid_area_cm2,first_pore_solid_area_cm2,first_three_pores_solid_area_cm2
```

Implemented:

- `precip_InitializePrecipitationAreaTimeseries.m`
- `precip_AppendPrecipitationAreaTimeseries.m`
- `precip_PNM_beauty3.m` initializes `precipitation_area_timeseries.csv`, appends one row per completed time step, and records the CSV path in `metadata.outputs`.
- `testPrecipitationAreaTimeseriesWritesRequiredColumns`
- `testPrecipitationAreaTimeseriesRejectsNonIntegerTimestep`

- [x] **Step 5: Export precipitation area plot**

Create:

```text
precipitation_area_timeseries.png
```

Plot:

```text
x-axis: time_min
y-axis: net solid area / precipitated area [cm2]
curves: total, first pore, first three pores
markers: 13, 18, 118 min
```

Implemented:

- `precip_PlotPrecipitationAreaTimeseries.m`
- `precip_PNM_beauty3.m` generates `precipitation_area_timeseries.png` after the run and records the PNG path in `metadata.outputs`.
- `testPrecipitationAreaTimeseriesPlotCreatesPng`

Smoke evidence:

```text
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260623_233235
```

This run generated:

```text
precipitation_area_timeseries.csv
precipitation_area_timeseries.png
```

This is export-wiring evidence only. The run ended at `finalTime_s = 0.015`, produced no 13/18/118 min benchmark snapshots, and still flagged `overshoot_c;mass_balance_drift`.

Review result:

- Parameter/source re-review passed for Step 4/5 because marker times are source-based but the short smoke is clearly framed as export wiring, not benchmark validation.
- Reviewer-risk re-review passed after hardening `precip_AppendPrecipitationAreaTimeseries.m` to reject non-integer timesteps and refreshing returned `result.metadata` after `metadata.final` is populated.

- [x] **Step 6: Export 13/18/118 min images**

Create directory:

```text
benchmark_snapshots/
```

For each target time, write:

```text
benchmark_snapshot_013min.png
benchmark_snapshot_018min.png
benchmark_snapshot_118min.png
```

Each image should show:

- Current solid/pore interface.
- Initial interface overlay.
- Split-inlet labels or equivalent annotation.
- Simulation time.

Implemented:

- `precip_BenchmarkSnapshotFilename.m` centralizes stable names:

```text
benchmark_snapshot_013min.png
benchmark_snapshot_018min.png
benchmark_snapshot_118min.png
benchmark_snapshot_001s.png
benchmark_snapshot_1p5s.png
```

- `precip_ExportBenchmarkSnapshots.m` exports PNG snapshots from triangular grid zero-level-set segments.
- `precip_PNM_beauty3.m` creates `benchmark_snapshots/`, records it as `metadata.outputs.benchmark_snapshots_dir`, calls the snapshot helper for due comparison times, and writes `metadata.benchmarkSnapshots`.
- Snapshot export failure is fatal for due benchmark targets; the run ending manifest errors with `RTSPHEM:Precipitate:MissingBenchmarkSnapshots` if a due nonempty PNG is missing.
- `precip_ConfigureZhangYoonBenchmark.m` records:

```matlab
cfg.benchmarkComparisonTimesSource = 'Zhang2010_Yoon2012_image_comparison_times';
cfg.geometryFraming = 'approximate_local_cylindrical_post_layout_not_digitized_geometry';
```

Tests added:

- `testBenchmarkSnapshotsSkipUntilDue`
- `testBenchmarkSnapshotsExportMinuteNamesWhenDue`, including `013min`, `018min`, and `118min`.
- `testBenchmarkSnapshotsUseSecondNamesForArtificialSmokeTimes`, including `001s`, `1p5s`, and `002s`.
- `testBenchmarkSnapshotsRejectInvalidGrid`
- `testBenchmarkSnapshotFilenameRejectsInvalidTarget`

- [x] **Step 7: Run export smoke test**

Run short benchmark with artificial comparison times:

```matlab
struct('endTime', 3, ...
       'timeStepperType', 'linear', ...
       'initialMacroscaleTimeStepSize', 1, ...
       'benchmarkComparisonTimes_s', [1 1.5 2], ...
       'benchmarkComparisonTimesSource', 'artificial_export_wiring_smoke', ...
       'numPartitionsMicroscale', 24, ...
       'meshNumPartitionsX', 24, ...
       'meshNumPartitionsY', 24, ...
       'reactionModel', 'tst', ...
       'enableNMRSimulation', false, ...
       'enableNMRSurrogate', false, ...
       'enablePNGSimulation', false)
```

Expected:

- `benchmark_comparison_times_log.csv` has three captured rows.
- `benchmark_snapshots` has three images.
- `precipitation_area_timeseries.csv` exists.
- `run_metadata.json` records `benchmarkComparisonTimesSource = artificial_export_wiring_smoke`.
- `metadata.benchmarkSnapshots.numTargets = 3`, `numExportedFiles = 3`, and `numMissingDue = 0`.

Observed smoke evidence:

```text
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260624_000244
```

Generated:

```text
benchmark_snapshots/benchmark_snapshot_001s.png
benchmark_snapshots/benchmark_snapshot_1p5s.png
benchmark_snapshots/benchmark_snapshot_002s.png
benchmark_comparison_times_log.csv
precipitation_area_timeseries.csv
precipitation_area_timeseries.png
run_metadata.json
```

This is artificial export-wiring evidence only. It does not validate Zhang/Yoon physics, does not produce real 13/18/118 min snapshots, and still showed transport overshoot / mass-balance diagnostics in the short run.

- [x] **Step 8: Review**

Parameter/source review: confirm 13/18/118 min are Zhang/Yoon image comparison times.  
Reviewer-risk review: confirm first-pore windows are marked as approximate unless using digitized geometry.

Review result:

- Parameter/source re-review passed: 13/18/118 min are recorded as Zhang/Yoon image comparison times, artificial smoke times are marked as `artificial_export_wiring_smoke`, and geometry framing is explicitly approximate / not digitized.
- Reviewer-risk re-review passed: snapshot export failures are no longer warning-only, run metadata includes a snapshot manifest, filenames are centralized, legend handles are stable, and tests cover 118 min plus decimal-second smoke naming.

Current status: comparison-time CSV, approximate area fields, reusable area helper, production area time-series CSV/PNG, stable benchmark snapshot helper, metadata snapshot manifest, and artificial snapshot smoke are complete. A later clipped diagnostic 13/18/118 min run produced real benchmark-time outputs and proved the export chain reaches the required comparison times, but physical Zhang/Yoon validation remains pending because that run used a finite PHREEQC transport limiter and retained CFL/mass-balance stability flags.

---

## Task 8: Compare with Zhang Experiment and Yoon Case 1/Case 5

**Files:**

- Create/modify: `ReactiveTransport/RTM/precipitate/precip_CompareZhangYoonBenchmark.m`
- Create/modify: `ReactiveTransport/RTM/precipitate/reference_data/zhang_yoon_reference_curves.csv`
- Create/modify: `ReactiveTransport/RTM/precipitate/reference_data/README.md`
- Test: `ReactiveTransport/RTM/precipitate/tests/test_precip_signed_helpers.m`
- Output: `zhang_yoon_area_comparison.csv`, `zhang_yoon_area_comparison.png`, `benchmark_comparison_report.md`

- [x] **Step 1: Prepare reference curve CSV schema**

Create:

```text
ReactiveTransport/RTM/precipitate/reference_data/zhang_yoon_reference_curves.csv
```

Columns:

```text
source,case,region,time_min,precipitated_area_norm,precipitated_area_cm2,note
```

Allowed `source` values:

```text
Zhang2010
Yoon2012
```

Allowed `case` values:

```text
experiment_25mM
case_1
case_5
```

Allowed `region` values:

```text
entire_domain
first_pore
first_three_pores
zhang_upgradient
zhang_middle
zhang_downgradient
```

Source-scoped region rules:

```text
Zhang2010 / experiment_25mM -> zhang_upgradient, zhang_middle, zhang_downgradient
Yoon2012 / case_1           -> entire_domain, first_pore, first_three_pores
Yoon2012 / case_5           -> entire_domain, first_pore, first_three_pores
```

Current repository state:

```text
source,case,region,time_min,precipitated_area_norm,precipitated_area_cm2,note
Zhang2010,experiment_25mM,zhang_upgradient,13,0.3833,NaN,digitized_approx|...
...
Yoon2012,case_5,first_pore,118,NaN,5.60e-5,digitized_approx|...
```

The file contains 18 approximate visual digitization rows: 9 Zhang selected-pore normalized pixel-area rows and 9 Yoon area rows converted from square micrometers to square centimeters. A temporary header-only CSV is still used in tests to verify that the comparison script rejects missing reference data.

- [x] **Step 1b: Add reference-data boundary README**

Create:

```text
ReactiveTransport/RTM/precipitate/reference_data/README.md
```

Required content:

```text
`zhang_yoon_reference_curves.csv` stores source-bounded approximate digitized reference rows.
See `digitization_notes.md` for axis calibration, unit conversion, and accuracy limits.
Do not add placeholder or synthetic values to this file.
```

- [x] **Step 2: Digitize initial reference data**

Use either existing extracted figure data or manual digitization from Zhang/Yoon plots.

Minimum required points:

```text
13 min
18 min
118 min
```

Preferred:

```text
full area-time curve points from the published figures
```

Record digitization note in the `note` column.

Current status: initial source-bounded visual digitization is complete for the required 13/18/118 min points. Local source assets exist under:

```text
literature_deng_reactive_transport/precipitationbench/extracted_figures/
```

Relevant candidate images include:

```text
zhang2010_si_page-4.png
yoon2012_fig4_area_si_ph_case1.png
yoon2012_fig7_case5_area_velocity_ph_si.png
```

Zhang 2010 source basis:

```text
Zhang Table 1: 25 mM is experiment 3 with Omega_c/Omega_v = 4.6/3.9.
Zhang Supporting Information Figure S3(b): selected pore space at Omega_c/Omega_v = 4.6/3.9.
Source-defined selected-pore regions: Upgradient, middle, Downgradient.
```

Do not use Zhang main-text Figure 5(b) as `experiment_25mM`; it is the 50 mM high-saturation experiment with Omega_c/Omega_v = 5.2/4.5.

Digitized outputs:

```text
ReactiveTransport/RTM/precipitate/reference_data/zhang_yoon_reference_curves.csv
ReactiveTransport/RTM/precipitate/reference_data/digitization_notes.md
```

Rows currently included:

```text
Zhang2010 / experiment_25mM / zhang_upgradient   -> 13, 18, 118 min
Zhang2010 / experiment_25mM / zhang_middle       -> 13, 18, 118 min
Zhang2010 / experiment_25mM / zhang_downgradient -> 13, 18, 118 min
Yoon2012 / case_1 / entire_domain                -> 13, 18, 118 min
Yoon2012 / case_5 / entire_domain                -> 13, 18, 118 min
Yoon2012 / case_5 / first_pore                   -> 13, 18, 118 min
```

Accuracy boundary:

```text
These are approximate visual digitizations from local rasterized figures.
They are suitable for first benchmark overlays and sanity checks.
They are not publisher tabular data and should be refined with WebPlotDigitizer or source tables before quantitative claims.
```

Do not enter synthetic rows into the repository reference CSV. Use temporary synthetic rows only inside unit tests.

- [x] **Step 3: Implement comparison script**

Function signature:

```matlab
function report = precip_CompareZhangYoonBenchmark(runDir, referenceCsv)
```

Inputs:

```text
runDir/precipitation_area_timeseries.csv
runDir/benchmark_comparison_times_log.csv
referenceCsv
```

Outputs:

```text
runDir/zhang_yoon_area_comparison.csv
runDir/zhang_yoon_area_comparison.png
runDir/benchmark_comparison_report.md
```

Implementation requirements now covered:

- Validate simulation input files exist:

```text
runDir/precipitation_area_timeseries.csv
runDir/benchmark_comparison_times_log.csv
```

- Validate reference CSV columns exactly include:

```text
source,case,region,time_min,precipitated_area_norm,precipitated_area_cm2,note
```

- Reject any reference CSV that has the correct schema but no digitized rows with `RTSPHEM:Precipitate:MissingReferenceCurves`. The repository reference CSV now contains source-bounded approximate digitized rows and is expected to pass this gate.

- Validate allowed source/case pairs:

```text
Zhang2010 -> experiment_25mM
Yoon2012  -> case_1
Yoon2012  -> case_5
```

- If reference rows have finite `precipitated_area_cm2` but missing `precipitated_area_norm`, fill normalized values per source/case/region using that group's maximum absolute cm2 value and document this in the report.

- State explicitly in the report that RTSPHEM normalized area is a simulation-internal normalization within each region by the simulation region's maximum absolute net area. It is not automatically the same normalization used in Zhang/Yoon figures.

- Include `entire_domain` because Yoon 2012 Case 1 and Case 5 report entire numerical-domain CaCO3 precipitate areas. Map RTSPHEM `total_net_solid_area_cm2` to the `entire_domain` simulation region when that column exists.

- Keep regions source-scoped. `Yoon2012 / case_1` and `Yoon2012 / case_5` may use `entire_domain`, `first_pore`, and `first_three_pores`. `Zhang2010 / experiment_25mM` may use only `zhang_upgradient`, `zhang_middle`, and `zhang_downgradient`, matching SI Figure S3(b)'s selected-pore regions.

- Preserve legacy run compatibility. If an older `precipitation_area_timeseries.csv` lacks `total_net_solid_area_cm2`, still compare `first_pore` and `first_three_pores`, omit RTSPHEM `entire_domain` rows, and write a report note.

- [x] **Step 4: Plot comparison scaffold**

Plot curves:

```text
RTSPHEM local benchmark
Zhang2010 experiment_25mM
Yoon2012 case_1
Yoon2012 case_5
```

Separate panels:

```text
entire domain
first pore
first three pores
Zhang upgradient
Zhang middle
Zhang downgradient
```

Current status: implemented and tested with both temporary synthetic reference rows and the repository source-bounded reference CSV. The script has also been run on the clipped diagnostic 13/18/118 min output directory, producing comparison CSV/PNG/report artifacts. Scientific Zhang/Yoon interpretation still requires a physically stable benchmark run, because the diagnostic run is intentionally flagged as output-chain evidence rather than quantitative validation.

- [x] **Step 5: Write report scaffold**

The report must include:

- Run directory and config.
- Geometry type: approximate array or digitized DXF/TIF.
- Flow and chemistry settings.
- Whether early precipitation growth is captured.
- Whether late area decrease or restructuring is captured.
- Which reference curve is closest.
- Remaining physical/modeling limitations.

Current implemented report scope:

- Writes `benchmark_comparison_report.md`.
- Reports the run directory and generated comparison files.
- Notes if cm2-only reference rows were normalized internally by group maximum.
- States that closest-reference and benchmark interpretation require higher-precision reference curves and a physically meaningful full run.

Full scientific interpretation is still pending because the digitization is approximate and the current 13/18/118 min comparison was generated from a clipped diagnostic run with unresolved transport stability flags.

- [x] **Step 6: Add comparison-script tests**

Tests added in:

```text
ReactiveTransport/RTM/precipitate/tests/test_precip_signed_helpers.m
```

Coverage:

```text
testCompareZhangYoonBenchmarkWritesOutputs
testCompareZhangYoonBenchmarkSupportsEntireDomainYoonCase1
testCompareZhangYoonBenchmarkRejectsZhangEntireDomainReference
testCompareZhangYoonBenchmarkSupportsZhangSelectedPoreRegions
testCompareZhangYoonBenchmarkRejectsYoonZhangSelectedPoreRegion
testCompareZhangYoonBenchmarkSupportsLegacyAreaCsvWithoutTotal
testCompareZhangYoonBenchmarkRejectsBadReferenceCsv
testCompareZhangYoonBenchmarkFillsReferenceNormFromCm2
testCompareZhangYoonBenchmarkRejectsInvalidSourceCasePair
testCompareZhangYoonBenchmarkAcceptsRepositoryReferenceCurves
testCompareZhangYoonBenchmarkRejectsHeaderOnlyReferenceCsv
testCompareZhangYoonBenchmarkReportsLimiterAndStabilityFlags
testComputePhreeqcTransportUpperBoundsUsesSplitInletMaxima
testComputePhreeqcTransportUpperBoundsAllowsInfiniteLimiter
testResolveAdaptiveTimeGridExtensionHonorsExplicitOverride
testResolveAdaptiveTimeGridExtensionUsesDefaultWhenUnset
testResolveAdaptiveTimeGridExtensionParsesFalseString
testResolveAdaptiveTimeGridExtensionRejectsNonScalar
testResolveAdaptiveTimeGridExtensionRejectsAmbiguousNumericAndText
testExtendTransportProblemsToStepperVisitsAllSpeciesAndVariables
testSignedBatchUsesPrecipLocalPhreeqcInputBuilder
testPrecipBuildCalcitePhreeqcInputSignedWritesSignedKinetics
```

Verification command:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate'); addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\couplePhreeqc'); results = runtests('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate\tests\test_precip_signed_helpers.m'); assertSuccess(results);"
```

Observed local result:

```text
60 tests passed
```

- [x] **Step 7: Static check comparison script**

Run:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "msg1 = checkcode('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate\precip_CompareZhangYoonBenchmark.m','-id'); msg2 = checkcode('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate\tests\test_precip_signed_helpers.m','-id'); assert(isempty(msg1)); assert(isempty(msg2));"
```

Observed local result: no `checkcode` messages.

- [x] **Step 8: Review**

Parameter/source review: verify reference data provenance and units.  
Reviewer-risk review: verify comparison does not overclaim strict reproduction when geometry is approximate.

Review result:

- Parameter/source review passed for the scaffold boundary when the repository CSV was still header-only. After source-bounded approximate digitized rows were added, a later parameter/source review passed for provenance, units, source/case/region scope, and the explicit approximate-digitization boundary.
- Reviewer-risk review initially failed because cm2-only reference rows were accepted without normalized plotting values, invalid source/case combinations could pass, and the report did not explain simulation-internal normalization.
- Fixes were applied: cm2-only rows are normalized per source/case/region group, source/case tuples are validated, the report explains RTSPHEM's internal normalization, and regression tests cover those cases.
- Reviewer-risk re-review passed. A later reviewer-risk pass confirmed the repository CSV smoke test now asserts 18 parsed reference rows, finite filled Yoon normalized values, scientific-notation cm2 parsing, Zhang selected-region presence, Zhang reference-only report messaging, and a temporary header-only rejection path.
- Later diagnostics/reporting review passed after the comparison report was updated to disclose finite `phreeqcTransportMaxFactor` use and observed `overshoot_c`, `advective_cfl_gt_1`, and `mass_balance_drift` flags.

Current status: comparison scaffold, source-bounded reference data, approximate digitization, diagnostic-run comparison outputs, and overclaim-prevention diagnostics are implemented. The repository reference CSV now contains approximate visual digitizations for 13/18/118 min, so it no longer raises `RTSPHEM:Precipitate:MissingReferenceCurves`. Quantitative benchmark interpretation still requires a physically stable run and preferably higher-precision WebPlotDigitizer/source-table reference data.

---

## Task 9: Full Benchmark Run

**Files:**

- Runner: `ReactiveTransport/RTM/precipitate/run_zhang_yoon_caco3_precipitation_benchmark.m`
- Config: `ReactiveTransport/RTM/precipitate/precip_ConfigureZhangYoonBenchmark.m`
- Outputs: `ReactiveTransport/RTM/precipitate/outputs/rtm_runs/<runName>/`

- [x] **Step 1: Run 25 mM local benchmark diagnostic**

Run:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate'); result = run_zhang_yoon_caco3_precipitation_benchmark(); disp(result.resultsDir);"
```

Expected:

- Run completes or stops with a documented numerical/PHREEQC error.
- Results are under `ReactiveTransport/RTM/precipitate/outputs/rtm_runs/<runName>/`.

Observed diagnostic run:

```text
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260624_095000
```

Run command summary:

```text
reactionModel = phreeqc
phreeqcRunGroup = phreeqc_database_calcite
phreeqcRateLaw = database_calcite
phreeqcDatabasePath = C:\Users\imgw\Downloads\RTSPHEM-P-main (1)\RTSPHEM-P-main\SourceCode\phreeqc-m.dat
endTime = 118 * 60 s
benchmarkComparisonTimes_s = [13 18 118] * 60
timeStepperType = linear
initialMacroscaleTimeStepSize = 60 s
meshNumPartitionsX = 24
meshNumPartitionsY = 24
enableNMRSimulation = false
enableNMRSurrogate = false
enablePNGSimulation = false
phreeqcTransportMaxFactor = 2
enableConcentrationCflLimit = false
```

This run completed to `finalTime_s = 7080` in `runtimeWallSeconds = 271.1650285`, but it is a clipped diagnostic/output-chain run, not a quantitative Zhang/Yoon validation.

Earlier unbounded preflight:

```text
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260624_012429
```

The unbounded preflight failed around `time_s = 60.62` with PHREEQC `A(H2O) Activity of water has not converged` after Ca/C/Na/Cl transport overshoot. This motivated `precip_ComputePhreeqcTransportUpperBounds.m`.

- [x] **Step 2: Inspect required outputs**

Required files:

```text
run_metadata.json
global_evolution_log.csv
phreeqc_summary_log.csv
benchmark_comparison_times_log.csv
precipitation_area_timeseries.csv
precipitation_area_timeseries.png
benchmark_snapshots/benchmark_snapshot_013min.png
benchmark_snapshots/benchmark_snapshot_018min.png
benchmark_snapshots/benchmark_snapshot_118min.png
```

Observed required outputs in `..._095000`:

```text
run_metadata.json
global_evolution_log.csv
phreeqc_results/phreeqc_summary_log.csv
benchmark_comparison_times_log.csv
precipitation_area_timeseries.csv
precipitation_area_timeseries.png
benchmark_snapshots/benchmark_snapshot_013min.png
benchmark_snapshots/benchmark_snapshot_018min.png
benchmark_snapshots/benchmark_snapshot_118min.png
```

Benchmark captures:

```text
13 min  -> timestep 13,  time_s = 780
18 min  -> timestep 18,  time_s = 1080
118 min -> timestep 118, time_s = 7080
```

- [ ] **Step 3: Validate physics sanity**

Check:

- pH remains finite.
- Ca/C/Na/Cl concentrations remain nonnegative or within expected numerical tolerance.
- `calcite_net_delta_moles_mean` has signed values.
- Precipitation area grows near mixing interface during early time.
- Permeability decreases when precipitation blocks pore space.

Observed diagnostic result:

```text
pH_mean remains finite around 6.91-6.95.
calcite_net_delta_moles_mean is signed but the mean is weakly negative at the sampled beginning/end.
benchmark comparison rows exist at 13/18/118 min.
permeability remained unchanged at 2.606564657559e9 mD.
porosity remained unchanged at 0.64575.
stability_diagnostics_log.csv reports overshoot_c;advective_cfl_gt_1;mass_balance_drift at the beginning and end.
advective_cfl is approximately 319.824 for the 60 s diagnostic timestep.
mass_balance_relative remains O(10^2).
```

Conclusion: physics sanity does not pass. This run is useful evidence that the export/comparison chain reaches 13/18/118 min, but it cannot be used for quantitative Zhang/Yoon validation.

- [x] **Step 4: Run comparison script**

Run:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate'); precip_CompareZhangYoonBenchmark('<runDir>', 'C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate\reference_data\zhang_yoon_reference_curves.csv');"
```

Expected:

```text
zhang_yoon_area_comparison.csv
zhang_yoon_area_comparison.png
benchmark_comparison_report.md
```

Observed outputs in `..._095000`:

```text
zhang_yoon_area_comparison.csv
zhang_yoon_area_comparison.png
benchmark_comparison_report.md
```

`benchmark_comparison_report.md` now contains a `Numerical Diagnostics` section:

```text
finite PHREEQC transport limiter used: phreeqcTransportMaxFactor = 2
Stability diagnostic flags observed: overshoot_c, advective_cfl_gt_1, mass_balance_drift
The run should be treated as diagnostic/output-chain evidence rather than quantitative Zhang/Yoon validation.
```

- [x] **Step 5: Review**

Parameter/source review: confirm final run uses Zhang/Yoon parameter set.  
Reviewer-risk review: inspect numerical stability, image overlays, area definitions, and comparison claims.

Review results:

- Parameter/source review passed for `precip_ComputePhreeqcTransportUpperBounds.m`: finite `phreeqcTransportMaxFactor` is scientifically defensible only as a numerical guard derived from configured initial and inlet chemistry; it is not a Zhang/Yoon literature parameter. Review caveat: any finite-cap run must report the factor and clipping boundary.
- Reviewer/code-risk review passed for the limiter: it addresses the immediate PHREEQC convergence failure by bounding H/Ca/C/Na/Cl before PHREEQC, but it clips mass and does not fix transport stability or mass balance.
- Parameter/source review passed for comparison-report diagnostics: the report correctly distinguishes literature/reference data from numerical limiter/stability diagnostics.
- Reviewer/code-risk review passed for comparison-report diagnostics: missing metadata/stability files are handled, finite limiter and stability flags are surfaced, and overclaiming is reduced. Caveats: malformed stability CSV is not caught, `ok` flags could be noisy, and future reports could quantify clipping counts.

- [ ] **Step 6: Produce physically stable benchmark run**

Required next run should avoid using clipping as the main stabilization mechanism. Candidate approaches:

```text
enableHardCflLimit = true
hardCflLimitInitialGrid = true or allowAdaptiveTimeGridExtension = true
initialMacroscaleTimeStepSize near the advective CFL candidate dt
coarser/faster exploratory mesh first, then refine
phreeqcTransportMaxFactor = Inf or sensitivity runs with factor values documented
```

Acceptance:

```text
No persistent advective_cfl_gt_1.
No persistent mass_balance_drift.
PHREEQC completes without convergence failure.
Precipitation area evolution is physically interpretable.
Permeability/porosity response is consistent with precipitation/clogging.
```

Current status: clipped diagnostic full-output run completed; physically stable benchmark remains pending.

Additional hard-CFL / adaptive-extension smoke evidence:

```text
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260624_101141
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260624_101830
```

Observed:

- `hardCflLimitInitialGrid = true` reduced the first-step advective CFL from the earlier 60 s diagnostic value of about 319.824 to 0.888.
- With `allowAdaptiveTimeGridExtension = false`, concentration/CFL shrinkage exhausted the preplanned grid and stopped around `finalTime_s = 0.3074`.
- With `allowAdaptiveTimeGridExtension = true`, the solver appended internal time steps and extended H/Ca/C/Na/Cl transport variables together, removing the previous "missing PHREEQC species data at appended step" failure.
- The no-clipping hard-CFL smoke still failed around `t = 0.336 s` because H concentration ran away and PHREEQC stopped with an `A(H2O) Activity of water has not converged` error. Therefore adaptive extension is a time-grid mechanics fix, not a physical stability fix.

Implemented after this smoke:

- `precip_ResolveAdaptiveTimeGridExtension.m` now strictly accepts scalar logical/numeric and true/false-like text, rejecting non-scalar or ambiguous values.
- `precip_ExtendTransportProblemsToStepper.m` provides a tested extension path for all active transport problems.
- `precip_PNM_beauty3.m` records adaptive append count, final time-grid summary, target end time, actual final time, and termination reason in metadata; PHREEQC exceptions now write `metadata.failure` before rethrowing.

---

## Task 10: Documentation and Handoff

**Files:**

- Modify: `literature_deng_reactive_transport/Zhang_Yoon_CaCO3_precipitation_benchmark_interface_plan.md`
- Create/modify: `literature_deng_reactive_transport/Zhang_Yoon_CaCO3_precipitation_execution_plan.md`
- Optional create: `ReactiveTransport/RTM/precipitate/README.md`

- [x] **Step 1: Document implementation boundary**

Document:

```text
All precipitation-specific code lives in ReactiveTransport/RTM/precipitate.
Original PNM_beauty3.m remains the dissolution solver.
```

- [x] **Step 2: Document completed code steps**

Record:

- Copied solver.
- Runner/config.
- Signed PHREEQC helpers.
- Signed interface-rate connection.
- Split inlet.
- Approximate Zhang/Yoon layout.
- Benchmark comparison-time CSV.
- Area time-series CSV/PNG export.
- Stable benchmark snapshot filenames and snapshot export.
- Comparison-script scaffold, source-scoped reference-data guard, and approximate digitized reference rows.
- PHREEQC transport upper-bound helper for finite numerical-guard runs.
- Diagnostic 13/18/118 min output-chain run and comparison-report stability disclosures.
- Existing tests and commands.

- [ ] **Step 3: Document full benchmark result**

After full run, add:

- Run path.
- MATLAB/COMSOL/PHREEQC environment.
- Whether PHREEQC COM was used.
- Number of time steps.
- Key output files.
- Comparison summary against Zhang/Yoon.

- [x] **Step 4: Add module README**

Create:

```text
ReactiveTransport/RTM/precipitate/README.md
```

Include:

```text
Purpose
How to run smoke tests
How to run local benchmark
How to interpret signed calcite delta
Known limitations
Output file descriptions
```

Current status: this execution plan is created and updated with completed code/scaffold steps through Task 9 diagnostic output-chain evidence, including approximate reference-curve digitization, `digitization_notes.md`, diagnostic comparison outputs, and the remaining physical-stability requirements. Module README, physically stable full-run result documentation, higher-precision reference refinement, and quantitative Zhang/Yoon comparison interpretation remain.
Current status: this execution plan is created and updated with completed code/scaffold steps through Task 9 diagnostic output-chain evidence, including approximate reference-curve digitization, `digitization_notes.md`, diagnostic comparison outputs, adaptive time-grid extension diagnostics, and the remaining physical-stability requirements. Module README exists. Physically stable full-run result documentation, higher-precision reference refinement, and quantitative Zhang/Yoon comparison interpretation remain.

---

## Review and Verification Policy

Every key implementation step must be followed by:

1. Parameter/source-basis review  
   Check whether values come from Zhang 2010, Yoon 2012, PHREEQC documentation, Deng-related precipitation papers, or are explicitly marked as implementation approximations.

2. Reviewer-perspective code-risk review  
   Check sign conventions, mass conservation, old-code isolation, numerical stability, output validity, and whether claims match what the code actually computes.

Minimum verification commands:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate'); addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\couplePhreeqc'); results = runtests('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate\tests\test_precip_signed_helpers.m'); assertSuccess(results);"
```

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "msg = checkcode('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate\precip_RunPhreeqcCalciteBatchSigned.m','-id'); assert(isempty(msg));"
```

Full benchmark verification requires PHREEQC COM and is not covered by helper tests.

Latest review checkpoint:

- Parameter/source subagent review passed for the split-inlet mesh-alignment safeguard and confirmed it is framed as a numerical implementation guard, not a Zhang/Yoon literature claim.
- Reviewer-risk subagent initially failed the mesh-alignment helper because tests covered only bare `N x 2` coordinates and not production-style `coordV` or `2 x N` coordinates.
- The requested fixes were applied: the helper now accepts `coordV` and explicit `2 x N` numeric coordinates, and tests cover both success and missing-endpoint failure paths.
- Parameter/source re-review passed after the fixes; no source-basis or smoke-framing issues remained.
- Reviewer-risk re-review passed for the narrow mesh-alignment/no-NMR smoke documentation change. Residual risks remain for benchmark physics, but they are not blockers for this guard/helper change.
- Task 7 Step 3 area-helper parameter/source review passed because the helper introduces no new literature parameter and explicitly marks first-pore and first-three-pore windows as approximate x-window diagnostics.
- Task 7 Step 3 reviewer-risk review passed after adding clipped-triangle, malformed-grid, and invalid level-set tests plus finite real grid/level validation. The helper is accepted as created/tested; production export wiring remains a later step.
- Task 7 Step 4/5 parameter/source review passed because the area time-series CSV/PNG are framed as local diagnostic exports; 13/18/118 min marker times are source-based, but the short smoke does not validate those comparison times.
- Task 7 Step 4/5 reviewer-risk review passed after adding integer timestep validation and refreshing returned `result.metadata` after `metadata.final`.
- Task 7 Step 6/7 parameter/source review initially failed because snapshot legends were ambiguous, artificial smoke times were not recorded in metadata, and the plan still described snapshot export as pending. Fixes added stable legend handles, `benchmarkComparisonTimesSource`, `geometryFraming`, metadata snapshot manifest, and updated plan status. Parameter/source re-review passed.
- Task 7 Step 6/7 reviewer-risk review initially failed because snapshot export failure was warning-only, metadata only recorded the snapshot directory, and tests did not cover 118 min or decimal-second filenames. Fixes made due snapshot failures fatal, added `metadata.benchmarkSnapshots`, centralized filename rules, and added 118 min / `1p5s` tests. Reviewer-risk re-review passed.
- Task 8 parameter/source review passed for the comparison scaffold and reference-data boundary when the repository CSV was header-only; later updates added source-bounded approximate digitized rows and are separately reviewed.
- Task 8 reviewer-risk review initially failed because cm2-only reference rows could pass without normalized plotting values, invalid source/case combinations were not rejected, and the report did not define simulation-internal normalization. Fixes added grouped cm2-to-normalized fill, source/case tuple validation, explicit normalization wording, and regression tests. Reviewer-risk re-review passed.
- Task 8 `entire_domain` source-boundary review initially failed because `entire_domain` is sourced from Yoon 2012 Fig. 4/Fig. 7 but was allowed for Zhang rows too, and the plan did not yet document the new region. Fixes made `entire_domain` valid only for `Yoon2012 / case_1` and `Yoon2012 / case_5`, added a Zhang `entire_domain` rejection test, documented the Yoon-only source basis, and updated the region/panel plan. Source-boundary re-review passed.
- Task 8 reviewer-risk review initially failed the `entire_domain` addition because requiring `total_net_solid_area_cm2` would break older area CSV outputs. Fixes made the total-area column optional, omitted RTSPHEM `entire_domain` rows for legacy CSVs, added a report note, and added a legacy compatibility test. Reviewer-risk re-review passed.
- Task 8 Zhang selected-pore schema update passed source-boundary review. The code allows `zhang_upgradient`, `zhang_middle`, and `zhang_downgradient` for `Zhang2010 / experiment_25mM`, rejects Yoon use of Zhang-specific regions, and records the source basis as Zhang Table 1 plus SI Figure S3(b).
- Task 8 Zhang selected-pore reviewer-risk review initially failed because generated reports did not explicitly say Zhang selected-pore panels are reference-only and have no matching RTSPHEM simulation rows. Fixes added report wording and a regression test asserting both phrases. Reviewer-risk re-review passed.
- Task 8 reference CSV digitization source/provenance review passed. The review confirmed Zhang 25 mM uses Table 1 experiment 3 and SI Figure S3(b), Yoon Case 1/Case 5 rows use Figure 4(a)/Figure 7(a), Zhang pixel units are kept normalized, Yoon square-micrometer units are converted to square centimeters, and all rows are marked as approximate visual digitizations. README wording was adjusted to distinguish Zhang normalized pixel rows from Yoon cm2 rows.
- Task 8 reference CSV reviewer-risk review initially failed because the header-only guard was no longer tested and the repository reference CSV test only checked for non-crashing behavior. Fixes restored a temporary header-only CSV rejection test and strengthened the repository reference CSV test to assert 18 parsed reference rows, finite filled Yoon normalized values, scientific-notation cm2 parsing, and Zhang reference-only report messaging. Reviewer-risk re-review passed.
- Task 9 PHREEQC transport-limiter source/parameter review passed. The review accepted `precip_ComputePhreeqcTransportUpperBounds.m` as a numerical guard derived from configured initial and split-inlet chemistry, not as a Zhang/Yoon literature parameter or calibration value.
- Task 9 PHREEQC transport-limiter reviewer-risk review passed. The review accepted the helper as an immediate PHREEQC convergence guard while noting that finite clipping can mask transport instability and cannot be used as proof of physical mass conservation.
- Task 9 comparison-report diagnostics source/parameter review passed. The report now separates source-bounded Zhang/Yoon reference curves from numerical limiter/stability diagnostics.
- Task 9 comparison-report diagnostics reviewer-risk review passed. The report now surfaces finite limiter use and `overshoot_c`, `advective_cfl_gt_1`, and `mass_balance_drift` flags, reducing the risk of overclaiming diagnostic output-chain results as quantitative validation.
- Task 9 adaptive time-grid extension source/parameter review passed for the first implementation: `allowAdaptiveTimeGridExtension` is a numerical runtime control, not a Zhang/Yoon parameter.
- Task 9 adaptive time-grid extension reviewer-risk review initially failed because the first implementation lacked direct species-extension test coverage, metadata for appended time grids, failure metadata, README wording for the no-clipping failure, and strict boolean parsing. Fixes added `precip_ExtendTransportProblemsToStepper.m`, stricter `precip_ResolveAdaptiveTimeGridExtension.m`, metadata/failure records, README wording, and tests. First re-review still failed because numeric values other than 0/1 and `yes/on/no/off` text were accepted, and the transport-extension test did not cover all variable fields. Fixes now restrict numeric/text parsing to `0/1/true/false` and test all 13 transport variable fields across five H/Ca/C/Na/Cl-like fixtures. Second reviewer-risk re-review passed.
- Task 9 signed PHREEQC input-builder isolation implemented by adding `precip_BuildCalcitePhreeqcInputSigned.m` and switching `precip_RunPhreeqcCalciteBatchSigned.m` to the precip-local builder. Source/parameter review passed because no new PHREEQC/literature parameter was introduced and the selected-output signed `KIN_DELTA("Calcite")` semantics were preserved. Reviewer-risk review initially failed because the isolation test was too weak, the test path setup could be shadowed, and copied error IDs still used the shared `RTSPHEM:Phreeqc:*` namespace. Fixes added stronger `which`/regexp tests, `addpath(moduleDir, '-begin')`, and precipitate-specific `RTSPHEM:Precipitate:*` error IDs. Reviewer-risk re-review passed.
- Task 9 PHREEQC output-helper isolation passed source/parameter review: the new precip-local CSV/plot helpers introduce no physical, Zhang/Yoon, or calibration parameter and only export existing signed PHREEQC diagnostics. Reviewer-risk review passed: copied PNM calls `precip_WritePhreeqcSpeciesTable` / `precip_ExportPhreeqcSpeciesPlots`, tests guard against bare shared-helper calls, signed CSV columns preserve parser/runner semantics, and local concentration plotting keeps the all-pore triangle mask convention.
- Task 9 split-inlet mass-balance diagnostic inlet term passed source/parameter review: the helper uses existing split inlet concentration fields, `yRange` lengths, velocity, and thickness; it introduces no Zhang/Yoon parameter or calibration value. Reviewer-risk review initially failed because the first fix risked leaving/removing the inlet concentration threshold ambiguously. The implementation now uses separate `diagnosticInletMassFlux` and `diagnosticInletConcentrationScale` helpers so mass-balance and overshoot/negative concentration flags retain distinct meanings. Reviewer-risk re-review passed.
- Task 9 concentration-CFL split-inlet scale unification passed source/parameter review: it introduces no new parameter and derives the overshoot scale from existing split-inlet boundary concentrations. Reviewer-risk review passed: the change only affects `ResolveConcentrationLimitedStep` time-step control, leaves transport/reaction/level-set paths unchanged, and is guarded by a regression test that prevents reverting the call to `initialHydrogenConcentration`.

---

## Remaining Gaps

1. Split-inlet helper tests and `SplitInletNotMeshAligned` guard tests are complete, including direct numeric coordinates, production-style `coordV`, missing-endpoint failure, and `2 x N` coordinate orientation. Parameter/source and reviewer-risk subagent re-reviews passed for this narrow change.
2. A no-NMR RTM-only smoke run completed with PNG-NMR, COMSOL NMR, and surrogate disabled. It remains wiring evidence only because it ended at `finalTime_s = 0.015` and still flagged `overshoot_c;mass_balance_drift`.
3. Full PNM signed-coupling is smoke verified and diagnostic-output verified. The chain reaches PHREEQC summaries, signed rates, benchmark-time snapshots, area-time CSV/PNG, and comparison outputs, but the physical precipitation/clogging response is not yet stable or benchmark-calibrated.
4. The PNG-enabled short run showed concentration overshoot, advective CFL greater than 1, and mass-balance drift. The clipped 13/18/118 min diagnostic run completed only after using `phreeqcTransportMaxFactor = 2` and disabling concentration/hard CFL stopping, and it still showed `overshoot_c;advective_cfl_gt_1;mass_balance_drift`. A physically valid benchmark needs tighter transport stabilization, concentration CFL control, adaptive time-grid handling, or another documented numerical treatment that removes these flags.
5. Area-time curve CSV/PNG and snapshot export are wired, short-smoke verified, and full diagnostic-output verified at 13/18/118 min. They should still be treated as output-chain evidence until produced by a physically stable run.
6. Stable 13/18/118 min snapshot filenames are implemented, tested, and generated in the clipped diagnostic run. Stable physical 13/18/118 min snapshots remain pending.
7. Zhang/Yoon reference CSV now contains approximate visual digitizations for required 13/18/118 min points; higher-precision WebPlotDigitizer/source-table refinement is still recommended before quantitative claims.
8. `zhang_yoon_area_comparison.csv`, `zhang_yoon_area_comparison.png`, and `benchmark_comparison_report.md` have been generated for the clipped diagnostic 13/18/118 min run. Equivalent outputs from a physically stable validation run are still pending.
9. Current Zhang/Yoon geometry is an approximate cylindrical-post layout; strict reproduction requires digitized DXF/TIF.
10. The worktree currently shows modifications to the original `ReactiveTransport/RTM/PNM_beauty3.m`; static audit found no direct call from the precipitation runner/solver to the original `PNM_beauty3.m`, but final delivery should still separate those original-solver edits from the precipitation-copy work.
11. The signed precipitation path now uses `precip_BuildCalcitePhreeqcInputSigned.m` instead of directly calling the modified shared `couplePhreeqc/BuildCalcitePhreeqcInput.m`. PHREEQC species CSV/plot output in the copied PNM now uses precip-local `precip_WritePhreeqcSpeciesTable.m` and `precip_ExportPhreeqcSpeciesPlots.m`. Source/parameter and reviewer-risk re-reviews passed for both isolation changes.

---

## Completed-Step Audit

This section records what has already been completed in the current working tree and what evidence supports it.

### Completed

- `ReactiveTransport/RTM/precipitate/precip_PNM_beauty3.m` exists as a copied solver with first-line function name `precip_PNM_beauty3(config)`.
- `ReactiveTransport/RTM/precipitate/run_zhang_yoon_caco3_precipitation_benchmark.m` calls the copied solver, not the original `PNM_beauty3.m`.
- `ReactiveTransport/RTM/precipitate/precip_ConfigureZhangYoonBenchmark.m` defines the local Zhang/Yoon benchmark configuration, including split CaCl2/Na2CO3 inlet fields and `mineralEvolutionMode = 'signed_calcite_surface'`.
- `ReactiveTransport/RTM/precipitate/precip_RunPhreeqcCalciteBatchSigned.m` and `precip_ParsePhreeqcSelectedOutputSigned.m` implement signed PHREEQC output handling.
- `ReactiveTransport/RTM/precipitate/precip_ScaleSignedCalciteDeltaToCellInventory.m`, `precip_ApplyPrescribedCalciteDissolutionSigned.m`, and `precip_ComputeSignedCalciteInterfaceRatePerArea.m` implement signed scaling/refill/rate conversion helpers.
- `ReactiveTransport/RTM/precipitate/precip_CreateTransportMultiInlet.m` implements the split-left inlet helper.
- `ReactiveTransport/RTM/precipitate/precip_RunLevelSetSignedVelocitySmoke.m` verifies the level-set sign convention against the real `LevelSetSolver2ndOrder` path.
- `ReactiveTransport/RTM/precipitate/tests/test_precip_signed_helpers.m` currently contains 69 MATLAB helper/integration tests, including split-inlet, mesh-alignment guard, precipitation area metrics, precipitation area CSV/plot export, benchmark snapshot export, comparison-script input/output guards, comparison-report limiter/stability diagnostics, PHREEQC transport upper-bound helper behavior, adaptive time-grid extension parsing, strict invalid boolean rejection, multi-species/all-variable transport extension, split-inlet diagnostic inlet flux/concentration-scale helpers, concentration-CFL split-inlet scale wiring, signed PHREEQC local-builder isolation, local PHREEQC input-builder error namespace checks, PHREEQC output-helper isolation, signed calcite CSV column export, local PHREEQC plot mask semantics, `entire_domain` Yoon Case 1 comparison support, Zhang `entire_domain` rejection, Zhang selected-pore region support, Yoon rejection of Zhang-specific regions, repository reference CSV semantic parsing, header-only reference rejection, legacy area-CSV compatibility, signed PHREEQC parser/scaler, PHREEQC COM directional cases, and level-set sign smoke coverage.
- `ReactiveTransport/RTM/precipitate/precip_BuildCalcitePhreeqcInputSigned.m` exists and is tested for signed PHREEQC input generation with `USER_PUNCH` and `KIN_DELTA("Calcite")`; `precip_RunPhreeqcCalciteBatchSigned.m` is tested to call this precip-local builder instead of the shared `BuildCalcitePhreeqcInput`.
- `ReactiveTransport/RTM/precipitate/precip_WritePhreeqcSpeciesTable.m`, `precip_ExportPhreeqcSpeciesPlots.m`, and `precip_PrepareConcentrationFaceData.m` exist and are tested. `precip_PNM_beauty3.m` is tested to call the precip-local PHREEQC output helpers instead of the shared `WritePhreeqcSpeciesTable` / `ExportPhreeqcSpeciesPlots` helpers.
- `ReactiveTransport/RTM/precipitate/precip_ComputePhreeqcTransportUpperBounds.m` exists and is tested. It computes optional PHREEQC transport upper bounds from initial, legacy inlet, and split-inlet chemistry fields and preserves the unbounded default when `phreeqcTransportMaxFactor = Inf`.
- `ReactiveTransport/RTM/precipitate/precip_ComputeBoundaryInletMassFlux.m` and `precip_ComputeBoundaryInletConcentrationScale.m` exist and are tested. Stability mass-balance diagnostics now use split-inlet segment fluxes for the inlet term while keeping a separate inlet concentration scale for `negative_c` / `overshoot_c` flags; `ResolveConcentrationLimitedStep` uses the same inlet concentration scale for concentration-CFL time-step control.
- `ReactiveTransport/RTM/precipitate/precip_ResolveAdaptiveTimeGridExtension.m` exists and is tested for default behavior, explicit override, false-string parsing, and non-scalar rejection.
- `ReactiveTransport/RTM/precipitate/precip_ExtendTransportProblemsToStepper.m` exists and is tested with five H/Ca/C/Na/Cl-like transport problem fixtures.
- `ReactiveTransport/RTM/precipitate/precip_ComputePrecipitationAreaMetrics.m` exists and is tested for total/window net area, Zhang-window defaults, clipped-triangle geometry, malformed grid rejection, and invalid level-set rejection.
- `ReactiveTransport/RTM/precipitate/precip_InitializePrecipitationAreaTimeseries.m`, `precip_AppendPrecipitationAreaTimeseries.m`, and `precip_PlotPrecipitationAreaTimeseries.m` exist and are tested for required CSV columns, integer timestep validation, and PNG creation.
- `ReactiveTransport/RTM/precipitate/precip_BenchmarkSnapshotFilename.m` and `precip_ExportBenchmarkSnapshots.m` exist and are tested for due/not-due export behavior, `013min` / `018min` / `118min` filenames, integer-second and decimal-second smoke filenames, invalid target validation, and malformed grid rejection.
- `ReactiveTransport/RTM/precipitate/precip_CompareZhangYoonBenchmark.m` exists and is tested with synthetic unit-test reference rows plus the repository reference CSV for output generation, `entire_domain` Yoon Case 1 support, Zhang `entire_domain` rejection, Zhang selected-pore support, legacy area-CSV compatibility, bad schema rejection, grouped cm2-to-normalized fill, invalid source/case/region rejection, and diagnostic reporting of finite PHREEQC transport limiters plus stability flags.
- `ReactiveTransport/RTM/precipitate/reference_data/zhang_yoon_reference_curves.csv` exists with 18 approximate visual digitization rows for 13/18/118 min.
- `ReactiveTransport/RTM/precipitate/reference_data/digitization_notes.md` documents source figures, approximate axis calibration, units, unit conversions, and accuracy limits.
- `ReactiveTransport/RTM/precipitate/reference_data/README.md` documents that no placeholder or synthetic reference values should be added and points readers to the digitization notes.
- `ReactiveTransport/RTM/precipitate/README.md` exists and documents module purpose, entry points, outputs, signed calcite interpretation, diagnostic-run boundaries, no-clipping hard-CFL failure evidence, and limitations.
- MATLAB helper/integration tests passed locally with:

```powershell
& "C:\Program Files\MATLAB\R2025b\bin\matlab.exe" -batch "addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate','-begin'); addpath('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\couplePhreeqc'); results = runtests('C:\Users\imgw\Documents\Codex\RTSPHEM-main\ReactiveTransport\RTM\precipitate\tests\test_precip_signed_helpers.m'); assertSuccess(results); fprintf('tests=%d passed=%d\n', numel(results), nnz([results.Passed]));"
```

- A short Zhang/Yoon local signed-PHREEQC smoke run completed and wrote:

```text
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260623_224132
```

- A no-NMR RTM-only short smoke run completed and wrote:

```text
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260623_225255
```

- A no-NMR area-export short smoke run completed and wrote:

```text
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260623_233235
```

It generated `precipitation_area_timeseries.csv` and `precipitation_area_timeseries.png`, but remains export-wiring evidence only because the run ended at `finalTime_s = 0.015` and still flagged `overshoot_c;mass_balance_drift`.

- A no-NMR artificial snapshot export smoke run completed and wrote:

```text
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260624_000244
```

It generated:

```text
benchmark_snapshots/benchmark_snapshot_001s.png
benchmark_snapshots/benchmark_snapshot_1p5s.png
benchmark_snapshots/benchmark_snapshot_002s.png
benchmark_comparison_times_log.csv
precipitation_area_timeseries.csv
precipitation_area_timeseries.png
run_metadata.json
```

The run metadata records `benchmarkComparisonTimesSource = artificial_export_wiring_smoke`, `geometryFraming = approximate_local_cylindrical_post_layout_not_digitized_geometry`, and `benchmarkSnapshots.numMissingDue = 0`. This is export-wiring evidence only; it is not a Zhang/Yoon physical validation.

- A clipped no-NMR 13/18/118 min diagnostic run completed and wrote:

```text
ReactiveTransport/RTM/precipitate/outputs/rtm_runs/zhang_yoon_caco3_25mM_20260624_095000
```

It generated:

```text
run_metadata.json
global_evolution_log.csv
phreeqc_results/phreeqc_summary_log.csv
benchmark_comparison_times_log.csv
precipitation_area_timeseries.csv
precipitation_area_timeseries.png
benchmark_snapshots/benchmark_snapshot_013min.png
benchmark_snapshots/benchmark_snapshot_018min.png
benchmark_snapshots/benchmark_snapshot_118min.png
zhang_yoon_area_comparison.csv
zhang_yoon_area_comparison.png
benchmark_comparison_report.md
```

The comparison-time log captured 13, 18, and 118 min exactly at `time_s = 780`, `1080`, and `7080`. The run is diagnostic output-chain evidence only because it used `phreeqcTransportMaxFactor = 2`, disabled concentration/hard CFL stopping, and retained `overshoot_c;advective_cfl_gt_1;mass_balance_drift` flags.

- A hard-CFL no-clipping short smoke without adaptive extension completed only to about `finalTime_s = 0.3074`, showing that CFL/concentration shrinkage could exhaust the preplanned time grid before `endTime`.
- A hard-CFL no-clipping short smoke with adaptive extension proved that appended steps can extend H/Ca/C/Na/Cl transport variables together, but it still failed around `t = 0.336 s` because H concentration ran away and PHREEQC reported `A(H2O) Activity of water has not converged`.

### Partly Completed

- Full PNM signed-coupling is partly verified: the short runs and clipped diagnostic run reached signed PHREEQC summaries and outputs, but the physical precipitation/clogging response is not yet stable or benchmark-calibrated.
- Benchmark comparison-time logging, area time-series plotting, and stable snapshot export are complete for a clipped diagnostic 13/18/118 min run. The same outputs still need to be regenerated from a physically stable benchmark run.
- `benchmark_comparison_times_log.csv` exists in both smoke outputs and the clipped diagnostic 13/18/118 min output. Only the latter reaches Zhang/Yoon comparison times, and it is explicitly diagnostic rather than physically validated.
- `precip_CompareZhangYoonBenchmark.m` is implemented as a comparison scaffold and data-boundary guard. It has been exercised with synthetic test rows, the repository approximate reference CSV, and a clipped diagnostic 13/18/118 min benchmark output.

### Not Completed

- Physically stable 13/18/118 min Zhang/Yoon validation run without relying on finite PHREEQC transport clipping as the primary stabilizer.
- Higher-precision Zhang experiment and Yoon Case 1/Case 5 reference curves beyond the current approximate visual digitization rows.
- Zhang/Yoon comparison outputs from a physically stable benchmark run.
- Final visual and curve comparison report for the stable benchmark.
- Final audit that original `ReactiveTransport/RTM/PNM_beauty3.m` changes are unrelated or intentionally separated.

---

## Completion Definition

This project is complete only when:

- `run_zhang_yoon_caco3_precipitation_benchmark.m` completes a local 25 mM benchmark run.
- The run produces 13/18/118 min images.
- The run produces precipitation area-time curves.
- The run is compared against Zhang experiment and Yoon Case 1/Case 5.
- Documentation states which outputs are strict comparisons and which are approximate because of geometry simplification.
- The precipitation module does not depend on unreviewed changes in original `ReactiveTransport/RTM/PNM_beauty3.m`; any original-solver modifications are either separated into their own task or explicitly documented as non-precipitation changes.

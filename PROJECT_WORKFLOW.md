# RTSPHEM Data-Generation Workflow

This repository is now organized around one main goal: generate RTM pore-scale
dissolution samples, then run COMSOL/NMR processing and T2 inversion for
NMR-agent training data.

## 1. Backup

A full backup was created before refactoring:

```text
C:\Users\imgw\Documents\Codex\RTSPHEM-main_backup_20260502_110923
```

## 2. Single RTM Run

MATLAB entry:

```text
ReactiveTransport\RTM\run_single_pnm.m
```

Edit the `cfg` struct near the top of that file:

- `cfg.layoutType`: `random`, `hex`, or `square`
- `cfg.characteristicLength`: random = target average throat, square/hex = minimum throat
- `cfg.inletVelocity`
- `cfg.initialHydrogenConcentration`
- output controls such as `cfg.exportEvery`, `cfg.saveIndividualPlots`, `cfg.saveFigureFiles`

The canonical solver is:

```text
ReactiveTransport\RTM\PNM_beauty3.m
```

It is now callable as:

```matlab
result = PNM_beauty3(cfg);
```

Default single-run outputs go to:

```text
outputs\rtm_runs\rtm_<timestamp>_<layout>\
```

Important files in each result folder:

- `run_metadata.json`: stable parameters and paths; do not rely on folder names for parameters
- `global_evolution_log.csv`
- `global_evolution.xlsx`
- `tortuosity_segments_log.csv`
- `tortuosity_segments.xlsx`
- `dxf_pore\`
- `dxf_solid\`
- `global_evolution_with_porosity_permeability.png`

## 3. Batch RTM Runs

MATLAB entry:

```text
ReactiveTransport\RTM\Batch_Simu\RunBatchExperiments.m
```

This still calls `PNM_batch`, but `PNM_batch` is now only a thin adapter around
`PNM_beauty3`. There is no longer a second solver implementation in the batch
path.

Default batch outputs go to:

```text
outputs\rtm_batches\batch_<timestamp>\exp_001\
outputs\rtm_batches\batch_<timestamp>\exp_002\
...
```

For NMR processing, keep:

```matlab
batchOutputOptions.exportEvery = 1;
batchOutputOptions.exportDXF = true;
```

For quick RTM debugging, you can reduce output cost:

```matlab
batchOutputOptions.exportEvery = 5;
batchOutputOptions.saveIndividualPlots = false;
batchOutputOptions.saveFigureFiles = false;
batchOutputOptions.saveRealtimePlot = false;
```

## 4. COMSOL/NMR Processing

Core COMSOL model:

```text
ReactiveTransport\NMR\CT-simulation.mph
```

Config file:

```text
ReactiveTransport\automation\AutomationConfig.m
```

Important fields:

- `data_root`: folder containing RTM result folders
- `mph_file`: COMSOL model file
- `python_exe`: Python interpreter for inversion
- `inversion_script`: `inversionBatch_matlab.py`
- `max_samples_per_folder`: number of DXF time steps to sample per RTM run
- `scale_factor`: DXF geometry scale factor used by COMSOL import

Batch NMR entry:

```text
ReactiveTransport\automation\run_automation.m
```

Single-folder NMR entry:

```text
ReactiveTransport\automation\run_single_folder.m
```

The automation now reads `run_metadata.json` first. Legacy
`dissolution_results-Da_...` folder names are only a fallback.

## 5. Cleanup Completed

The following paths were removed after explicit confirmation. They are still
available in the full backup listed above.

During MATLAB smoke testing, two grid classes were found to be required by the
RTM solver and restored from the backup into `ReactiveTransport\src\LevelSetSolver2ndOrder`:

```text
CartesianGrid.m
FoldedCartesianGrid.m
```

```text
3DPermCNN
3DRelDiffCNN
Data
LevelSetSolverOld
LevelSetTriplePoint
ReactiveTransport\Machine Learning
ReactiveTransport\littleScripts
ReactiveTransport\scripts
ReactiveTransport\automation\__pycache__
ReactiveTransport\HyPHM\doc
ReactiveTransport\RTM\BatchConvertdxf.m
ReactiveTransport\RTM\BatchConvertdxf_v2.m
ReactiveTransport\automation\quick_test.m
ReactiveTransport\automation\test_long_path.m
ReactiveTransport\automation\README_长路径修复.md
ReactiveTransport\automation\长路径处理修复说明.md
```

Keep these core paths:

```text
ReactiveTransport\RTM
ReactiveTransport\NMR
ReactiveTransport\automation
ReactiveTransport\src
ReactiveTransport\HyPHM\classes
ReactiveTransport\HyPHM\tools
ReactiveTransport\HyPHM\domains
ReactiveTransport\HyPHM\opt
ReactiveTransport\HyPHM\symbolic
ReactiveTransport\ExternalRoutines
```

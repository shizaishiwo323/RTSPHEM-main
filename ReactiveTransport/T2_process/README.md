# nmr_t2 standard package

This directory has been standardized as a Python package for NMR T2 processing.

## Features

- **NNLS inversion** with MATLAB-compatible core behavior
- **L-curve regularization selection** with reciprocal-slope criterion
- **Paired plotting** for raw decay and T2 spectrum
- **Gaussian peak decomposition** for spectrum interpretation

## Package structure

- `nmr_t2/config.py`: typed configuration objects
- `nmr_t2/io_utils.py`: robust Excel parsing and standardized exports
- `nmr_t2/nnls.py`: fixed-regularization inversion
- `nmr_t2/lcurve.py`: L-curve inversion
- `nmr_t2/plotting.py`: figure generation
- `nmr_t2/gaussian.py`: Gaussian decomposition
- `nmr_t2/pipelines.py`: high-level workflows

## RTSPHEM automation usage

The RTSPHEM automation bridge uses fixed-regularization NNLS for numerical
simulation outputs:

- regularization / smoothing factor: `0.01`
- no L-curve search is used in the automation path
- visualization is generated through `nmr_t2.plotting`

## Standard output naming

All high-level pipelines export files under this convention:

`<dataset_name>__<artifact_name>.<ext>`

Examples:

- `SimulationDecay__nnls_spectrum.xlsx`
- `SimulationDecay__nnls_summary.csv`
- `SimulationDecay__lcurve_metrics.xlsx`
- `SimulationDecay__gaussian_summary.xlsx`

## Quick start

You can run the complete workflow demo from `main.ipynb`.

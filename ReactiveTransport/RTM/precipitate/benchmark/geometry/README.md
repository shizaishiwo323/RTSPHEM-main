# Zhang/Yoon Geometry Packages

Place calibrated Zhang/Yoon geometry packages here. Each package should have a
`geometry_manifest.json` file validated by
`precip_LoadZhangYoonGeometryPackage`.

Required package assets:

- source micromodel image used for geometry extraction
- axis or length calibration CSV
- processing script that creates the masks
- substrate mask file (`.mat` or row/column `.csv`)
- first-pore and first-three-pore region mask file (`.mat` or row/column `.csv`)
- uncertainty CSV for calibration and mask extraction

Required manifest fields:

```json
{
  "packageName": "zhang_yoon_geometry",
  "sourceFigure": "Zhang/Yoon micromodel image",
  "sourceImageFile": "source_micromodel.png",
  "calibrationCsv": "axis_calibration.csv",
  "processingScript": "build_masks.m",
  "substrateMaskFile": "substrate_mask.mat",
  "regionMaskFile": "region_masks.mat",
  "uncertaintyCsv": "geometry_uncertainty.csv"
}
```

Do not commit placeholder geometry images or guessed masks. Until a package
passes the validator, the benchmark should keep reporting approximate geometry
and non-quantitative status.

The validator records `assetFilesVerified = true` only after all required
assets exist and both substrate/region masks parse successfully. Readiness also
requires positive substrate and first-pore mask cell counts, with the
first-three-pores mask at least as large as the first-pore mask.

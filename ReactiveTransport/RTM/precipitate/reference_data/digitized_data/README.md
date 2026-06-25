# Digitized Reference Data Packages

Place source-locked Zhang/Yoon digitization packages here. Each package should
have its own subdirectory and a `digitization_manifest.json` file validated by
`precip_LoadReferenceDigitizationPackage`.

Required package assets:

- source screenshot used for digitization
- WebPlotDigitizer project file
- axis-calibration CSV
- raw digitizer export CSV
- conversion script from raw export to repository reference schema
- uncertainty CSV
- converted reference CSV with the standard columns

Do not commit placeholder screenshots, synthetic digitizer projects, or guessed
reference values. Until a package passes the validator, reference curves remain
non-quantitative.

# Bundled PHREEQC Databases for RTM Coupling

This directory keeps the PHREEQC database files used by the RTM PHREEQC
coupling scripts so benchmark runs do not depend on files under a user
Downloads directory.

- `database/phreeqc_rates.dat`: preferred database for calcite
  kinetic/speciation coupling in the RTSPHEM PHREEQC path. This file uses
  the official PHREEQC `Calcite` rate law with `PARM(1)` in `cm^2/mol
  calcite` and `PARM(2)` as the `M/M0` surface-area exponent.
- `database/phreeqc.dat`: fallback database from the same official PHREEQC
  database snapshot.
- `database/backup_legacy_20260626/`: previous local `phreeqc-m.dat` and
  `phreeqc.dat` copies retained for reproducibility.

Only the database files are bundled here. The MATLAB implementation still
calls `IPhreeqcCOM.Object`, so the USGS IPhreeqcCOM engine must remain
installed and registered on Windows.

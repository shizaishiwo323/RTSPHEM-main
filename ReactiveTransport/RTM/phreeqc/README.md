# Bundled PHREEQC Databases for RTM Coupling

This directory keeps the PHREEQC database files used by the RTM PHREEQC
coupling scripts so benchmark runs do not depend on files under a user
Downloads directory.

- `database/phreeqc-m.dat`: preferred database for calcite kinetic/speciation
  coupling in the RTSPHEM PHREEQC path.
- `database/phreeqc.dat`: fallback database copied from the same reference
  source.

Only the database files are bundled here. The MATLAB implementation still
calls `IPhreeqcCOM.Object`, so the USGS IPhreeqcCOM engine must remain
installed and registered on Windows.

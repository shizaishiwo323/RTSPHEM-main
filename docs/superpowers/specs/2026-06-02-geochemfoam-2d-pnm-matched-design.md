# GeoChemFoam 2D PNM-Matched Simulation Design

## Objective

Create a new GeoChemFoam-based workflow for a two-dimensional microfluidic calcite dissolution simulation that matches the geometry intent and input parameters of `ReactiveTransport/RTM/run_single_pnm.m`.

The workflow must preserve the existing PNM data and outputs, place all new files under a dedicated directory, run through Docker when the Docker Linux engine is available, and record enough metadata for later comparison with PNM/NMR-agent dataset preparation.

## New Directory

All new simulation assets will be placed under the user-requested external root:

```text
C:\Users\imgw\Documents\Codex\geochemfoam
```

Planned subdirectories:

```text
C:\Users\imgw\Documents\Codex\geochemfoam\
  README.md
  docker/
  scripts/
  geometry/
  cases/
    pnm_matched_2d_calcite_ale/
  manifests/
  logs/
```

This keeps the GeoChemFoam workflow separate from `ReactiveTransport/RTM/`, `ReactiveTransport/NMR/`, existing `outputs/rtm_runs/`, and `Data/`.

## Source Template Choice

GeoChemFoam provides a clear two-dimensional tutorial:

```text
tutorials/reactiveTransport/reactiveTransportFoam/2DGrainstoneBimolecular
```

This case uses `n_z=1` and `frontandback` patches with `empty`, which matches the required OpenFOAM treatment for a 2D simulation.

GeoChemFoam also provides a calcite dissolution and moving-interface tutorial:

```text
tutorials/reactiveTransport/reactiveTransportALEFoam/3DcalcitePostALE
```

This case uses `reactiveTransportALEFoam`, `reactiveSurfaceConcentrationMixed`, and `pointMotionU` to move dissolving calcite surfaces. It is physically closer to `run_single_pnm.m` than the 2D bimolecular tutorial.

The new case will combine these ideas:

- Use the 2D mesh and `empty` boundary structure from `2DGrainstoneBimolecular`.
- Use the calcite dissolution fields and ALE surface motion from `3DcalcitePostALE`.
- Keep `reactiveTransportALEFoam` as the primary solver for the first matched workflow.

## Docker Route

The workflow will use Docker as the execution environment. The local machine currently has Docker installed, but the Docker Desktop Linux engine was not running during discovery. The implementation will therefore include:

- A Docker build script.
- A Docker run wrapper that mounts `C:\Users\imgw\Documents\Codex\geochemfoam` into the container.
- Environment checks that report Docker availability before compilation or simulation.
- Logs for build, solver compilation, mesh creation, and case execution.

The workflow will not assume a WSL Ubuntu distribution is installed.

## Geometry Input

The PNM script currently uses:

```text
C:\Users\imgw\Documents\Codex\论文复现\pore-scale-simulation-reproduction\organized_gpu_ac3d_reproduction_20260527\data\dissolution_results-Da_40.4424_Pe_4.1640_L_0.1200_square\validation-domin.dxf
```

PNM DXF settings:

- Domain layers: `domin`, `DOMAIN`
- Solid layer: `calcite`
- Reference length: `1500` DXF units equals `0.15 cm`
- Import direction: `rotate90_cw`
- Current raster resolution: `200 x 100`

The GeoChemFoam workflow will convert the same DXF geometry into a 2D OpenFOAM-compatible solid boundary input. The first implementation will generate an extruded thin STL from the calcite solid mask, because GeoChemFoam tutorials already use STL plus `snappyHexMesh` for image-derived pore geometries.

The front/back thickness will be a small numerical thickness only. The OpenFOAM case will remain 2D by using one cell in the empty direction and `empty` front/back patches.

## Parameter Mapping

The workflow will read or duplicate the current PNM parameters and write an explicit manifest. Unit conversions will be written in the generated configuration rather than hidden in scripts.

| Quantity | PNM value | GeoChemFoam/OpenFOAM value |
| --- | ---: | ---: |
| Inlet velocity | `0.1 cm/s` | `1.0e-3 m/s` |
| Diffusion coefficient | `5.0e-5 cm^2/s` | `5.0e-9 m^2/s` |
| Inlet H+ concentration | `1.37e-5 mol/cm^3` | `13.7 mol/m^3` |
| Characteristic length | `0.001 cm` | `1.0e-5 m` |
| Domain length X | `0.2662331517 cm` from metadata | `2.662331517e-3 m` |
| Domain length Y | `0.15 cm` from metadata | `1.5e-3 m` |
| End time | `16200 s` | `16200 s` |

The PNM metadata from the recent matched run records:

- `Da = 0.00101106`
- `Pe = 2`
- `flowDirection = left_to_right`
- `layoutType = external_dxf`

Reaction-rate handling must be explicit because `run_single_pnm.m` and GeoChemFoam ALE boundary fields use different parameter conventions. The implementation will create a documented `reaction_mapping.json` and write the value used for GeoChemFoam `k` in the manifest. If exact dimensional equivalence cannot be proven from the solver dictionaries, the manifest will mark this as a calibrated/assumed mapping rather than exact.

## Case Flow

The generated workflow will follow this order:

1. Check Docker status.
2. Clone or update `shizaishiwo323/GeoChemFoam-5.2` inside `ReactiveTransport/GeoChemFoam2D/GeoChemFoam-5.2`.
3. Build or prepare the Docker image.
4. Compile GeoChemFoam solvers inside the container.
5. Convert the PNM DXF geometry into a 2D extruded STL and record geometry metadata.
6. Create `cases/pnm_matched_2d_calcite_ale/` from GeoChemFoam tutorial assets.
7. Patch case dictionaries for the matched domain, mesh, inlet/outlet, empty front/back boundaries, transport properties, concentration, reaction boundary, and run time.
8. Run mesh generation.
9. Run the flow/reactive transport solver.
10. Extract comparable outputs into `manifests/` and `logs/`.

## Outputs

Expected outputs:

- Case directory with OpenFOAM dictionaries and initial fields.
- `geometry/pnm_matched_2d_calcite_ale.stl`
- `geometry/pnm_matched_geometry.json`
- `manifests/pnm_matched_2d_calcite_ale_manifest.json`
- `logs/docker_build.log`
- `logs/geochemfoam_compile.log`
- `logs/create_mesh.log`
- `logs/reactive_transport.log`
- `outputs/process_images/velocity/velocity_t*.png`
- `outputs/process_images/concentration/concentration_t*.png`
- `outputs/process_images/interface/interface_t*.png`
- `outputs/process_images/combined/combined_t*.png`
- OpenFOAM time directories in the case folder if the solver runs successfully.

The implementation will not overwrite PNM `.mat`, `.mph`, DXF, Excel, COMSOL, or inversion outputs.

## Validation

Before claiming the workflow is complete, the implementation must check:

- Docker is installed and whether the Linux engine is running.
- The GeoChemFoam repository is available at the expected commit or branch.
- The source DXF exists.
- Geometry conversion produced a non-empty STL.
- The case has `frontandback` patches with `empty`.
- Mesh generation produced `constant/polyMesh`.
- If solver time directories exist, the process-image script can generate velocity, concentration, interface, and combined process PNGs similar in purpose to `outputs/rtm_runs/rtm_20260602_184459_192_external_dxf/timestep_0001.png`.
- The case manifest records PNM parameters, converted OpenFOAM parameters, geometry source, script version, creation date, and success/failure status.

If Docker Desktop is not running, the implementation can still generate the directory, scripts, case templates, geometry converter, and manifests, but must clearly report that solver compilation and simulation were not run.

## Scope Boundaries

This design does not modify `ReactiveTransport/RTM/run_single_pnm.m`.

This design does not alter existing PNM outputs under `outputs/rtm_runs/`.

This design does not run COMSOL, T2 inversion, or NMR surrogate processing.

This design focuses on a single matched 2D GeoChemFoam case first. Batch dataset generation can be added only after the single case is verified.

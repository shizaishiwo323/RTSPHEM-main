# GeoChemFoam 2D PNM-Matched Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an external Docker-oriented GeoChemFoam 2D calcite dissolution workflow rooted at `C:\Users\imgw\Documents\Codex\geochemfoam`.

**Architecture:** The external root contains configuration, Python utilities, PowerShell wrappers, Docker files, generated geometry, case assets, manifests, logs, and process images. Python scripts handle deterministic tasks that can be verified locally: PNM parameter mapping, DXF-to-thin-STL conversion, case dictionary generation, and process-image generation. Docker wrappers handle GeoChemFoam clone/build/run when Docker Desktop's Linux engine is available.

**Tech Stack:** Windows PowerShell, Python standard library `unittest`, Docker, OpenFOAM/GeoChemFoam shell scripts, GeoChemFoam-5.2 from `shizaishiwo323/GeoChemFoam-5.2`.

---

### Task 1: Scaffold External Workflow Root

**Files:**
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\README.md`
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\config\pnm_matched_2d_calcite_ale.json`
- Create directories: `docker`, `scripts`, `tests`, `geometry`, `cases`, `manifests`, `logs`, `outputs\process_images\velocity`, `outputs\process_images\concentration`, `outputs\process_images\interface`, `outputs\process_images\combined`

- [ ] **Step 1: Create directories**

Run:

```powershell
$root = 'C:\Users\imgw\Documents\Codex\geochemfoam'
$dirs = @(
  $root,
  "$root\config",
  "$root\docker",
  "$root\scripts",
  "$root\tests",
  "$root\geometry",
  "$root\cases",
  "$root\manifests",
  "$root\logs",
  "$root\outputs\process_images\velocity",
  "$root\outputs\process_images\concentration",
  "$root\outputs\process_images\interface",
  "$root\outputs\process_images\combined"
)
foreach ($dir in $dirs) { New-Item -ItemType Directory -Force -Path $dir | Out-Null }
```

Expected: directories exist, no existing data is deleted.

- [ ] **Step 2: Write root README**

Create `README.md` explaining the Docker route, 2D case, directory structure, and command sequence:

```text
powershell -ExecutionPolicy Bypass -File .\scripts\check_environment.ps1
powershell -ExecutionPolicy Bypass -File .\scripts\clone_geochemfoam.ps1
python .\scripts\setup_case.py --config .\config\pnm_matched_2d_calcite_ale.json
powershell -ExecutionPolicy Bypass -File .\scripts\build_docker.ps1
powershell -ExecutionPolicy Bypass -File .\scripts\run_in_docker.ps1 -Command compile
powershell -ExecutionPolicy Bypass -File .\scripts\run_in_docker.ps1 -Command mesh
powershell -ExecutionPolicy Bypass -File .\scripts\run_in_docker.ps1 -Command run
python .\scripts\generate_process_images.py --config .\config\pnm_matched_2d_calcite_ale.json
```

- [ ] **Step 3: Write config JSON**

Create `config\pnm_matched_2d_calcite_ale.json` with these values:

```json
{
  "case_name": "pnm_matched_2d_calcite_ale",
  "root": "C:\\Users\\imgw\\Documents\\Codex\\geochemfoam",
  "geochemfoam_repo": "https://github.com/shizaishiwo323/GeoChemFoam-5.2.git",
  "geochemfoam_branch": "main",
  "source_dxf": "C:\\Users\\imgw\\Documents\\Codex\\论文复现\\pore-scale-simulation-reproduction\\organized_gpu_ac3d_reproduction_20260527\\data\\dissolution_results-Da_40.4424_Pe_4.1640_L_0.1200_square\\validation-domin.dxf",
  "domain_layers": ["domin", "DOMAIN"],
  "solid_layers": ["calcite"],
  "dxf_reference_length_units": 1500.0,
  "dxf_reference_length_cm": 0.15,
  "dxf_import_direction": "rotate90_cw",
  "domain_length_x_cm": 0.26623315171331874,
  "domain_length_y_cm": 0.15,
  "empty_thickness_m": 1.0e-6,
  "mesh_nx": 600,
  "mesh_ny": 338,
  "mesh_nz": 1,
  "pnm_parameters": {
    "Da": 0.00101106,
    "Pe": 2.0,
    "inlet_velocity_cm_s": 0.1,
    "diffusion_cm2_s": 5.0e-5,
    "inlet_h_mol_cm3": 1.37e-5,
    "characteristic_length_cm": 0.001,
    "rate_coefficient_tst_mol_dm2_s": 1.0e-4,
    "molar_volume_cm3_mol": 36.9,
    "end_time_s": 16200,
    "write_interval_s": 750
  },
  "geochemfoam_parameters": {
    "inlet_velocity_m_s": 0.001,
    "diffusion_m2_s": 5.0e-9,
    "inlet_h_mol_m3": 13.7,
    "domain_length_x_m": 0.0026623315171331874,
    "domain_length_y_m": 0.0015,
    "reaction_k_assumed": 1.0e-7,
    "stoich_coeff": 1.0,
    "solid_density_kg_m3": 2710.0,
    "solid_mw_kg_kmol": 100.0,
    "kinematic_viscosity_m2_s": 1.0e-6,
    "end_time_s": 16200,
    "write_interval_s": 750,
    "initial_delta_t_s": 0.01,
    "max_delta_t_s": 1.0
  }
}
```

- [ ] **Step 4: Verify scaffold**

Run:

```powershell
Test-Path 'C:\Users\imgw\Documents\Codex\geochemfoam\config\pnm_matched_2d_calcite_ale.json'
```

Expected: `True`.

### Task 2: Add Tests for Config and Geometry Conversion

**Files:**
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\tests\test_workflow.py`
- Create later: `C:\Users\imgw\Documents\Codex\geochemfoam\scripts\workflow_common.py`
- Create later: `C:\Users\imgw\Documents\Codex\geochemfoam\scripts\convert_dxf_to_stl.py`

- [ ] **Step 1: Write failing tests**

Create `tests\test_workflow.py` with tests that import missing modules:

```python
import json
import tempfile
import unittest
from pathlib import Path

from scripts.workflow_common import cm_to_m, load_config, paths_for_case
from scripts.convert_dxf_to_stl import DxfPolyline, write_extruded_stl


ROOT = Path(r"C:\Users\imgw\Documents\Codex\geochemfoam")


class WorkflowConfigTests(unittest.TestCase):
    def test_cm_to_m_converts_values(self):
        self.assertAlmostEqual(cm_to_m(0.15), 0.0015)

    def test_config_contains_matched_pnm_values(self):
        cfg = load_config(ROOT / "config" / "pnm_matched_2d_calcite_ale.json")
        self.assertEqual(cfg["case_name"], "pnm_matched_2d_calcite_ale")
        self.assertAlmostEqual(cfg["geochemfoam_parameters"]["inlet_velocity_m_s"], 0.001)
        self.assertAlmostEqual(cfg["geochemfoam_parameters"]["diffusion_m2_s"], 5.0e-9)
        self.assertAlmostEqual(cfg["geochemfoam_parameters"]["inlet_h_mol_m3"], 13.7)

    def test_paths_for_case_use_external_root(self):
        cfg = load_config(ROOT / "config" / "pnm_matched_2d_calcite_ale.json")
        paths = paths_for_case(cfg)
        self.assertEqual(paths["root"], ROOT)
        self.assertEqual(paths["case_dir"], ROOT / "cases" / "pnm_matched_2d_calcite_ale")


class GeometryTests(unittest.TestCase):
    def test_write_extruded_stl_creates_ascii_solid(self):
        polyline = DxfPolyline(
            layer="calcite",
            closed=True,
            points=[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        )
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "solid.stl"
            write_extruded_stl([polyline], output, scale=0.001, thickness=1e-6)
            text = output.read_text(encoding="utf-8")
        self.assertIn("solid pnm_matched_calcite", text)
        self.assertIn("facet normal", text)
        self.assertIn("endsolid pnm_matched_calcite", text)


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run tests to verify RED**

Run:

```powershell
cd C:\Users\imgw\Documents\Codex\geochemfoam
python -m unittest tests.test_workflow -v
```

Expected: FAIL with `ModuleNotFoundError` for `scripts.workflow_common`.

### Task 3: Implement Shared Workflow Utilities and DXF-to-STL Conversion

**Files:**
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\scripts\__init__.py`
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\scripts\workflow_common.py`
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\scripts\convert_dxf_to_stl.py`

- [ ] **Step 1: Implement `workflow_common.py`**

Functions:

```python
def cm_to_m(value):
    return float(value) * 0.01

def load_config(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))

def paths_for_case(config):
    root = Path(config["root"])
    case_name = config["case_name"]
    return {
        "root": root,
        "repo": root / "GeoChemFoam-5.2",
        "case_dir": root / "cases" / case_name,
        "geometry_dir": root / "geometry",
        "manifest_dir": root / "manifests",
        "logs_dir": root / "logs",
        "outputs_dir": root / "outputs",
    }

def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
```

- [ ] **Step 2: Implement `convert_dxf_to_stl.py`**

Implement a minimal ASCII DXF LWPOLYLINE parser for group codes `8`, `10`, `20`, `70`; transform `rotate90_cw`; scale DXF units to meters; write a thin ASCII STL with top, bottom, and side triangles for solid-layer polylines.

- [ ] **Step 3: Run tests to verify GREEN**

Run:

```powershell
cd C:\Users\imgw\Documents\Codex\geochemfoam
python -m unittest tests.test_workflow -v
```

Expected: all tests pass.

### Task 4: Add Docker and Windows Wrapper Scripts

**Files:**
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\docker\Dockerfile`
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\docker\container_entrypoint.sh`
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\scripts\check_environment.ps1`
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\scripts\clone_geochemfoam.ps1`
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\scripts\build_docker.ps1`
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\scripts\run_in_docker.ps1`

- [ ] **Step 1: Write Dockerfile**

Use `opencfd/openfoam-default:2212` when available, install Python tooling needed by GeoChemFoam tutorials, and set `/workspace` as the mounted root.

- [ ] **Step 2: Write `container_entrypoint.sh`**

Support commands:

```bash
compile
mesh
run
postprocess
shell
```

Each command writes logs under `/workspace/logs`.

- [ ] **Step 3: Write PowerShell wrappers**

The wrappers check Docker status, clone the repo, build the image, and run a chosen command. They must not delete directories.

- [ ] **Step 4: Verify wrapper syntax**

Run:

```powershell
powershell -ExecutionPolicy Bypass -File C:\Users\imgw\Documents\Codex\geochemfoam\scripts\check_environment.ps1
```

Expected: reports Docker client and Docker engine status without modifying data.

### Task 5: Generate Case Files and Manifests

**Files:**
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\scripts\setup_case.py`
- Generated: `C:\Users\imgw\Documents\Codex\geochemfoam\geometry\pnm_matched_2d_calcite_ale.stl`
- Generated: `C:\Users\imgw\Documents\Codex\geochemfoam\geometry\pnm_matched_geometry.json`
- Generated: `C:\Users\imgw\Documents\Codex\geochemfoam\cases\pnm_matched_2d_calcite_ale\...`
- Generated: `C:\Users\imgw\Documents\Codex\geochemfoam\manifests\pnm_matched_2d_calcite_ale_manifest.json`

- [ ] **Step 1: Write `setup_case.py`**

This script loads config, converts DXF to STL, writes `system/blockMeshDict`, `system/controlDict`, `constant/transportProperties`, `constant/thermoPhysicalProperties`, initial fields `0/C`, `0/U`, `0/p`, `0/pointMotionU`, and a manifest.

- [ ] **Step 2: Run setup**

Run:

```powershell
cd C:\Users\imgw\Documents\Codex\geochemfoam
python .\scripts\setup_case.py --config .\config\pnm_matched_2d_calcite_ale.json
```

Expected: generated STL, case dictionaries, and manifest exist.

- [ ] **Step 3: Verify generated case**

Run:

```powershell
Test-Path C:\Users\imgw\Documents\Codex\geochemfoam\geometry\pnm_matched_2d_calcite_ale.stl
Select-String -Path C:\Users\imgw\Documents\Codex\geochemfoam\cases\pnm_matched_2d_calcite_ale\system\blockMeshDict -Pattern 'frontandback|empty'
Test-Path C:\Users\imgw\Documents\Codex\geochemfoam\manifests\pnm_matched_2d_calcite_ale_manifest.json
```

Expected: first and third return `True`; second shows both `frontandback` and `empty`.

### Task 6: Add Process Image Generation

**Files:**
- Create: `C:\Users\imgw\Documents\Codex\geochemfoam\scripts\generate_process_images.py`
- Generated if solver outputs exist: `C:\Users\imgw\Documents\Codex\geochemfoam\outputs\process_images\...\*.png`

- [ ] **Step 1: Write `generate_process_images.py`**

Use Python standard library plus `tkinter` when available to create PNG summaries from OpenFOAM scalar/vector field files when solver outputs exist. If solver time directories are absent, the script writes a clear status into the manifest and exits with code `0` without claiming images were generated.

- [ ] **Step 2: Run image script**

Run:

```powershell
cd C:\Users\imgw\Documents\Codex\geochemfoam
python .\scripts\generate_process_images.py --config .\config\pnm_matched_2d_calcite_ale.json
```

Expected: if no solver time directories exist, manifest is updated with `process_images.status = "no_solver_outputs"`.

### Task 7: Try Docker-Aware Execution

**Files:**
- No new files unless Docker is available.
- Logs under: `C:\Users\imgw\Documents\Codex\geochemfoam\logs`

- [ ] **Step 1: Run environment check**

Run:

```powershell
powershell -ExecutionPolicy Bypass -File C:\Users\imgw\Documents\Codex\geochemfoam\scripts\check_environment.ps1
```

Expected: reports whether Docker engine is available.

- [ ] **Step 2: Clone GeoChemFoam repository**

Run:

```powershell
powershell -ExecutionPolicy Bypass -File C:\Users\imgw\Documents\Codex\geochemfoam\scripts\clone_geochemfoam.ps1
```

Expected: repository exists or is reported as already present.

- [ ] **Step 3: Build and run only if Docker engine is available**

Run:

```powershell
powershell -ExecutionPolicy Bypass -File C:\Users\imgw\Documents\Codex\geochemfoam\scripts\build_docker.ps1
powershell -ExecutionPolicy Bypass -File C:\Users\imgw\Documents\Codex\geochemfoam\scripts\run_in_docker.ps1 -Command compile
powershell -ExecutionPolicy Bypass -File C:\Users\imgw\Documents\Codex\geochemfoam\scripts\run_in_docker.ps1 -Command mesh
powershell -ExecutionPolicy Bypass -File C:\Users\imgw\Documents\Codex\geochemfoam\scripts\run_in_docker.ps1 -Command run
```

Expected: if Docker engine is unavailable, wrappers stop with a clear message and logs explain the blocked step. If Docker engine is available, logs show progress or the first GeoChemFoam/OpenFOAM error.

### Task 8: Final Verification

**Files:**
- Read generated files and logs only.

- [ ] **Step 1: Run local tests**

Run:

```powershell
cd C:\Users\imgw\Documents\Codex\geochemfoam
python -m unittest tests.test_workflow -v
```

Expected: all tests pass.

- [ ] **Step 2: Verify required artifacts**

Run:

```powershell
$root = 'C:\Users\imgw\Documents\Codex\geochemfoam'
Test-Path "$root\README.md"
Test-Path "$root\geometry\pnm_matched_2d_calcite_ale.stl"
Test-Path "$root\cases\pnm_matched_2d_calcite_ale\system\blockMeshDict"
Test-Path "$root\manifests\pnm_matched_2d_calcite_ale_manifest.json"
```

Expected: all return `True`.

- [ ] **Step 3: Report execution status**

Final report must state which data-preparation stage this affects: GeoChemFoam 2D RTM alternative, DXF geometry conversion, Docker environment setup, solver execution attempt, and process-image generation. It must also state whether MATLAB, COMSOL, Python inversion, and GeoChemFoam solver runs actually executed.

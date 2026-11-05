# dipole_calculus

**`dipole_calculus`** is a modular Python package for calculating spatial differential properties — **Curl** ($\nabla \times \vec{p}$), **Divergence** ($\nabla \cdot \vec{p}$), **Gradient** ($\nabla \|\vec{p}\|$), and **Laplacian** ($\nabla^2 \vec{p}$) — of water electric dipoles across GROMACS molecular dynamics trajectories. It supports multi-core parallel execution, vector/scalar analysis modes, Gaussian `.cube` volumetric exports, and automatic VMD `.vmd` visualization script generation.

---

## Project Structure

```text
dipole_calculus/
├── pyproject.toml                  # Package configuration and metadata (PEP 621)
├── requirements.txt                # Runtime dependencies
├── README.md                       # Project documentation
├── main.py                         # Direct CLI runner script
└── src/
    └── dipole_calculus/
        ├── __init__.py             # Public package API
        ├── cli.py                  # Command-line interface & result formatting
        ├── core/
        │   ├── __init__.py
        │   ├── calculator.py       # Main trajectory & multiprocessing orchestrator
        │   ├── grid.py             # 3D meshgrid construction & differential operators
        │   └── worker.py           # Parallel chunk worker for ProcessPoolExecutor
        └── io/
            ├── __init__.py
            ├── exporters.py        # Gaussian .cube and VMD .vmd script exporters
            └── topology.py         # Recursive GROMACS .top / .itp charge parser
```

---

## Differential Properties (`--property`)

Given the interpolated 3D dipole vector field $\vec{V}(x, y, z) = (V_x, V_y, V_z)$:

1. **`curl`** ($\nabla \times \vec{V}$) — *3D Vector Field*:
   - Measures local rotational vorticity of the electric dipole field.
   - Supports `--mode magnitude`, `--mode vector`, and `--mode dot_product`.
2. **`divergence`** ($\nabla \cdot \vec{V}$) — *3D Scalar Field*:
   - Measures local dipole flux sources ($>0$) and sinks ($<0$).
   - Inherently a scalar field; if `--mode` is passed, a non-blocking warning is displayed and standard scalar calculation proceeds.
3. **`gradient`** ($\nabla \|\vec{V}\|$) — *3D Vector Field*:
   - Computes the spatial gradient of the dipole magnitude $V_{\text{mag}} = \sqrt{V_x^2 + V_y^2 + V_z^2}$.
   - Supports `--mode magnitude`, `--mode vector`, and `--mode dot_product`.
4. **`laplacian`** ($\nabla^2 \vec{V}$) — *3D Vector Field*:
   - Computes the component-wise vector Laplacian $(\nabla^2 V_x, \nabla^2 V_y, \nabla^2 V_z)$.
   - Supports `--mode magnitude`, `--mode vector`, and `--mode dot_product`.

---

## Installation

### 1. Editable installation (Recommended for development)
```bash
pip install -e .
```
Once installed, the `dipole-calculus` command is available directly in your terminal.

### 2. Install dependencies only
```bash
pip install -r requirements.txt
```

---

## Usage

You can run the tool either via the `main.py` entry-point script or via the installed `dipole-calculus` command:

```bash
python main.py \
    -top topol.top \
    -ref system.gro \
    -traj trajectory.xtc \
    -box 40.0 40.0 40.0 \
    -g 80 80 80 \
    -cpus 8 \
    -box_ref "protein" \
    -sel "resname SOL" \
    --property curl \
    --mode magnitude \
    --gen_tcl
```

### Output Files
Output files are dynamically named according to the selected `--property`:
- **Time-averaged volumetric grid**: `mean_{property}.cube` (e.g., `mean_curl.cube`, `mean_divergence.cube`)
- **Per-frame volumetric grid** (when `--cube_per_frame N` is set): `frame_{frame_idx}_{property}.cube`
- **VMD visualization script** (when `--gen_tcl` is set): `plot_{property}.vmd`

### Command-Line Arguments

| Flag | Description | Default |
| :--- | :--- | :--- |
| `-top`, `--topol` | Path to the global GROMACS topology file (`.top`) | *Required* |
| `-ref`, `--reference` | Path to the reference coordinate file (`.gro`) | *Required* |
| `-traj`, `--trajectory` | Path to the trajectory file (`.xtc`) | *Required* |
| `-box`, `--box` | Total grid box dimensions along X, Y, Z (in Å) | *Required* |
| `-g`, `--grid` | Number of grid points along X, Y, Z | `80 80 80` |
| `-cpus`, `--parallel` | Number of CPU workers for parallel execution | `os.cpu_count()` |
| `-box_ref`, `--box_ref` | MDAnalysis selection string for the grid center | `"protein"` |
| `-sel`, `--sel` | MDAnalysis selection string for target solvent molecules | `"resname SOL"` |
| `--property` | Differential operator: `curl`, `divergence`, `gradient`, or `laplacian` | `"curl"` |
| `--mode` | Calculation mode for vector operators: `magnitude`, `vector`, or `dot_product` | `"magnitude"` |
| `--cube_per_frame` | Frame interval to export individual `frame_N_{property}.cube` files | `None` |
| `--gen_tcl` | Generate `plot_{property}.vmd` script for VMD visualization | `False` |

---

## License

This project is dual-licensed under both the **MIT License** and the **GNU General Public License v3.0 or later (GPL-3.0-or-later)**. You may choose either license at your option:

- For the permissive terms, see [LICENSE-MIT](LICENSE-MIT).
- For the copyleft terms, see [LICENSE-GPL](LICENSE-GPL).
- For details on the dual-licensing structure, see [LICENSE](LICENSE).


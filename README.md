# Dipole calculus

**`dipole_calculus`** is a modular Python package for calculating spatial differential properties — **Curl** ($\nabla \times \vec{p}$), **Divergence** ($\nabla \cdot \vec{p}$), **Gradient** ($\nabla \|\vec{p}\|$), and **Laplacian** ($\nabla^2 \vec{p}$) — of solvent (Ex.: water) electric dipoles across GROMACS molecular dynamics trajectories. It supports multi-core parallel execution, vector/scalar analysis modes, Gaussian `.cube` volumetric exports, and automatic VMD `.vmd` visualization script generation.

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

## Analysis Modes (`--mode`)

For vector properties (**`curl`**, **`gradient`**, and **`laplacian`**), the `--mode` flag determines how the 3D vector $\vec{A} = (A_x, A_y, A_z)$ is processed, reported, and exported:

### 1. `magnitude` (Default)
- **Mathematical Definition**: 
$\|\vec{A}\| = \sqrt{A_x^2 + A_y^2 + A_z^2}$
- **Physical Meaning**: Quantifies the absolute intensity of the field at every grid voxel, regardless of its spatial orientation.
- **Output**: 
  - The terminal prints the spatial mean magnitude for each trajectory frame.
  - The generated `.cube` file contains a 3D scalar volume representing the temporal average of the field magnitude across the simulation box (ideal for isosurface visualization in VMD).

### 2. `vector`
- **Mathematical Definition**: Preserves the full three-dimensional vector $\vec{A} = (A_x, A_y, A_z)$ at every grid point ($3 \times N_x \times N_y \times N_z$).
- **Physical Meaning**: Retains complete directional information and orientation. It allows you to verify whether the field possesses a preferential net polarization or alignment along a specific Cartesian axis ($X$, $Y$, or $Z$) across time.
- **Output**:
  - The terminal outputs individual per-axis averages: $\langle A_x \rangle$, $\langle A_y \rangle$, $\langle A_z \rangle$ for each frame.
  - The generated `.cube` file saves the magnitude of the temporally averaged vector field $\sqrt{\langle A_x \rangle^2 + \langle A_y \rangle^2 + \langle A_z \rangle^2}$.

  

### 3. `dot_product` (Radial Projection)
- **Mathematical Definition**: 
$$\vec{A} \cdot \hat{r} \quad \text{where} \quad \hat{r} = \frac{\vec{r} - \vec{r}_{\text{center}}}{\|\vec{r} - \vec{r}_{\text{center}}\|}$$
is the unit position vector pointing outward from the reference center of mass (e.g., protein center).
- **Physical Meaning**: Measures the projection of the vector field along the radial direction relative to the solute:
  - **Positive values ($> 0$)**: The vector points outward (away from the reference center).
  - **Negative values ($< 0$)**: The vector points inward (toward the reference center).
  - Highly useful for investigating whether dipole vorticity or gradient fluxes are radially oriented around a macromolecule.
- **Output**:
  - The terminal prints the frame-averaged scalar dot product.
  - The generated `.cube` file exports the 3D signed scalar field showing inward vs. outward orientation.

> [!NOTE]
> **Summary on Divergence vs Modes**: Since **divergence** is already an intrinsic scalar quantity ($\nabla \cdot \vec{V} \in \mathbb{R}$), the `--mode` flag does not apply to it. Divergence is always processed and exported as a signed scalar field where positive regions indicate dipole flux expansion and negative regions indicate convergence.


## Installation

### 1. Editable installation
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

```bash
dipole-calculus \
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


# Dipole calculus

**`dipole_calculus`** is a modular Python package for calculating spatial differential properties — **Curl** (∇ × p), **Divergence** (∇ · p), **Gradient** (∇|p|), and **Laplacian** (∇²p) — of solvent (e.g., water) electric dipoles across GROMACS molecular dynamics trajectories. It supports multi-core parallel execution, vector/scalar analysis modes, Gaussian `.cube` volumetric exports, and automatic VMD `.vmd` visualization script generation.

---

## Project Structure

```text
dipole_calculus/
├── pyproject.toml                  # Package configuration and metadata (PEP 621)
├── requirements.txt                # Runtime dependencies
├── README.md                       # Project documentation
├── main.py                         # Direct CLI runner script
├── dipole_calculus_standalone.py   # Single-file standalone script (HPC-ready)
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

Given the interpolated 3D dipole vector field V(x, y, z) = (Vx, Vy, Vz):

1. **`curl`** — *3D Vector Field*:
   - Measures local rotational vorticity of the electric dipole field:

$$
\text{curl}_x = \frac{\partial V_z}{\partial y} - \frac{\partial V_y}{\partial z}, \quad
\text{curl}_y = \frac{\partial V_x}{\partial z} - \frac{\partial V_z}{\partial x}, \quad
\text{curl}_z = \frac{\partial V_y}{\partial x} - \frac{\partial V_x}{\partial y}
$$

   - Supports `--mode magnitude`, `--mode vector`, and `--mode dot_product`.

2. **`divergence`** — *3D Scalar Field*:
   - Measures net volumetric flux divergence: local sources of dipole expansion (> 0) and sinks/convergence (< 0):

$$
\text{div} = \nabla \cdot \vec{V} = \frac{\partial V_x}{\partial x} + \frac{\partial V_y}{\partial y} + \frac{\partial V_z}{\partial z}
$$

   - **Mode behavior**: Divergence is inherently a scalar quantity and has no directional components (X, Y, Z) to decompose or project. Therefore, `--mode` is **not needed**. If `--mode` is passed, a non-blocking terminal warning is displayed and standard scalar calculation proceeds automatically.

3. **`gradient`** — *3D Vector Field*:
   - Computes the spatial variation vector of the local dipole magnitude:

$$
\text{grad}_i = \frac{\partial V_{\text{mag}}}{\partial x_i} \quad \text{where} \quad V_{\text{mag}} = \sqrt{V_x^2 + V_y^2 + V_z^2}
$$

   - Supports `--mode magnitude`, `--mode vector`, and `--mode dot_product`.

4. **`laplacian`** — *3D Vector Field*:
   - Computes the component-wise spatial diffusion/curvature vector:

$$
\nabla^2 V_i = \frac{\partial^2 V_i}{\partial x^2} + \frac{\partial^2 V_i}{\partial y^2} + \frac{\partial^2 V_i}{\partial z^2} \quad \text{for } i \in \{x, y, z\}
$$

   - Supports `--mode magnitude`, `--mode vector`, and `--mode dot_product`.

---

## Analysis Modes (`--mode`)

For vector properties (**`curl`**, **`gradient`**, and **`laplacian`**), the `--mode` flag determines how the 3D vector A = (Ax, Ay, Az) is processed, reported, and exported:

### 1. `magnitude` (Default)
- **Mathematical Definition**:

$$
\|\vec{A}\| = \sqrt{A_x^2 + A_y^2 + A_z^2}
$$

- **Physical Meaning**: Quantifies the absolute intensity of the field at every grid voxel, regardless of its spatial orientation.
- **Output**: 
  - The terminal prints the spatial mean magnitude for each trajectory frame.
  - The generated `.cube` file contains a 3D scalar volume representing the temporal average of the field magnitude across the simulation box (ideal for isosurface visualization in VMD).

### 2. `vector`
- **Mathematical Definition**: Preserves the complete 3D vector A = (Ax, Ay, Az) across every grid point (3 x Nx x Ny x Nz).
- **Physical Meaning**: Retains complete directional information and orientation. It allows you to verify whether the field possesses a preferential net polarization or alignment along a specific Cartesian axis (X, Y, or Z) across time.
- **Output**:
  - The terminal outputs individual per-axis averages: ⟨Ax⟩, ⟨Ay⟩, ⟨Az⟩ for each frame.
  - The generated `.cube` file saves the magnitude of the temporally averaged vector field:

$$
\|\langle \vec{A} \rangle\| = \sqrt{\langle A_x \rangle^2 + \langle A_y \rangle^2 + \langle A_z \rangle^2}
$$

### 3. `dot_product` (Radial Projection)
- **Mathematical Definition**: Radial scalar projection of the vector field A onto the radial unit vector r̂:

$$
\vec{A} \cdot \hat{r} \quad \text{where} \quad \hat{r} = \frac{\vec{r} - \vec{r}_c}{\| \vec{r} - \vec{r}_c \|}
$$

  Here, rc denotes the reference center of mass.
- **Physical Meaning**: Measures the projection of the vector field along the radial direction relative to the solute:
  - **Positive values (> 0)**: The vector points outward (away from the reference center).
  - **Negative values (< 0)**: The vector points inward (toward the reference center).
  - Highly useful for investigating whether dipole vorticity or gradient fluxes are radially oriented around a macromolecule.
- **Output**:
  - The terminal prints the frame-averaged scalar dot product.
  - The generated `.cube` file exports the 3D signed scalar field showing inward vs. outward orientation.

> [!NOTE]
> **Summary on Divergence vs Modes**: Since **divergence** is already an intrinsic scalar quantity (∇ · V), the `--mode` flag does not apply to it. Divergence is always processed and exported as a signed scalar field where positive regions indicate dipole flux expansion and negative regions indicate convergence.

---

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

"""Spatial grid generation and electric dipole differential operator routines."""

import numpy as np
from scipy.interpolate import griddata


def build_spatial_grid(box_center, box_edges, grid_resolution):
    """Build 3D meshgrid, bounding coordinates, and step sizes around the box center."""
    radius_x = box_edges[0] / 2.0
    radius_y = box_edges[1] / 2.0
    radius_z = box_edges[2] / 2.0

    min_x, max_x = box_center[0] - radius_x, box_center[0] + radius_x
    min_y, max_y = box_center[1] - radius_y, box_center[1] + radius_y
    min_z, max_z = box_center[2] - radius_z, box_center[2] + radius_z

    min_coords = [min_x, min_y, min_z]
    max_coords = [max_x, max_y, max_z]

    x_coords = np.linspace(min_x, max_x, grid_resolution[0])
    y_coords = np.linspace(min_y, max_y, grid_resolution[1])
    z_coords = np.linspace(min_z, max_z, grid_resolution[2])

    grid_x, grid_y, grid_z = np.meshgrid(x_coords, y_coords, z_coords, indexing="ij")

    dx = x_coords[1] - x_coords[0]
    dy = y_coords[1] - y_coords[0]
    dz = z_coords[1] - z_coords[0]

    return grid_x, grid_y, grid_z, dx, dy, dz, min_coords, max_coords


def _compute_scalar_laplacian(scalar_grid, dx, dy, dz):
    """Compute the 3D scalar Laplacian via two successive gradient passes."""
    df_dx, df_dy, df_dz = np.gradient(scalar_grid, dx, dy, dz)
    d2f_dx2 = np.gradient(df_dx, dx, axis=0)
    d2f_dy2 = np.gradient(df_dy, dy, axis=1)
    d2f_dz2 = np.gradient(df_dz, dz, axis=2)
    return d2f_dx2 + d2f_dy2 + d2f_dz2


def _apply_vector_mode(vec_x, vec_y, vec_z, grid_x, grid_y, grid_z, box_center, mode):
    """Reduce or project a 3D vector field according to magnitude, vector, or dot_product mode."""
    if mode == "vector":
        frame_mean = (float(np.mean(vec_x)), float(np.mean(vec_y)), float(np.mean(vec_z)))
        frame_contribution = np.stack([vec_x, vec_y, vec_z], axis=0)
        export_grid = np.sqrt(vec_x**2 + vec_y**2 + vec_z**2)
    elif mode == "dot_product":
        # Compute position vectors relative to the center of mass
        x_rel = grid_x - box_center[0]
        y_rel = grid_y - box_center[1]
        z_rel = grid_z - box_center[2]

        # Normalize radial position vectors (avoiding division by zero)
        r_mag = np.sqrt(x_rel**2 + y_rel**2 + z_rel**2)
        r_mag[r_mag == 0] = 1e-10

        u_x = x_rel / r_mag
        u_y = y_rel / r_mag
        u_z = z_rel / r_mag

        # Dot product between vector field and unit radial position vector
        dot_grid = (vec_x * u_x) + (vec_y * u_y) + (vec_z * u_z)

        frame_mean = float(np.mean(dot_grid))
        frame_contribution = dot_grid
        export_grid = dot_grid
    else:
        mag_grid = np.sqrt(vec_x**2 + vec_y**2 + vec_z**2)
        frame_mean = float(np.mean(mag_grid))
        frame_contribution = mag_grid
        export_grid = mag_grid

    return frame_mean, frame_contribution, export_grid


def compute_frame_property(
    oxygen_positions,
    dipole_vectors,
    grid_x,
    grid_y,
    grid_z,
    dx,
    dy,
    dz,
    box_center,
    property_name="curl",
    mode="magnitude",
):
    """Interpolate dipoles onto the 3D grid and compute the selected differential property."""
    # Interpolate dipole components onto the 3D grid using nearest-neighbor method
    vx_grid = griddata(oxygen_positions, dipole_vectors[:, 0], (grid_x, grid_y, grid_z), method="nearest")
    vy_grid = griddata(oxygen_positions, dipole_vectors[:, 1], (grid_x, grid_y, grid_z), method="nearest")
    vz_grid = griddata(oxygen_positions, dipole_vectors[:, 2], (grid_x, grid_y, grid_z), method="nearest")

    if property_name == "divergence":
        # Compute spatial divergence: dVx/dx + dVy/dy + dVz/dz
        dvx_dx = np.gradient(vx_grid, dx, axis=0)
        dvy_dy = np.gradient(vy_grid, dy, axis=1)
        dvz_dz = np.gradient(vz_grid, dz, axis=2)
        div_grid = dvx_dx + dvy_dy + dvz_dz

        frame_mean = float(np.mean(div_grid))
        return frame_mean, div_grid, div_grid

    if property_name == "gradient":
        # Compute spatial gradient of the dipole vector magnitude: grad(|V|)
        v_mag = np.sqrt(vx_grid**2 + vy_grid**2 + vz_grid**2)
        grad_x, grad_y, grad_z = np.gradient(v_mag, dx, dy, dz)
        return _apply_vector_mode(grad_x, grad_y, grad_z, grid_x, grid_y, grid_z, box_center, mode)

    if property_name == "laplacian":
        # Compute component-wise vector Laplacian: (lap(Vx), lap(Vy), lap(Vz))
        lap_x = _compute_scalar_laplacian(vx_grid, dx, dy, dz)
        lap_y = _compute_scalar_laplacian(vy_grid, dx, dy, dz)
        lap_z = _compute_scalar_laplacian(vz_grid, dx, dy, dz)
        return _apply_vector_mode(lap_x, lap_y, lap_z, grid_x, grid_y, grid_z, box_center, mode)

    # Default to curl: nabla x V
    dvx_dx, dvx_dy, dvx_dz = np.gradient(vx_grid, dx, dy, dz)
    dvy_dx, dvy_dy, dvy_dz = np.gradient(vy_grid, dx, dy, dz)
    dvz_dx, dvz_dy, dvz_dz = np.gradient(vz_grid, dx, dy, dz)

    curl_x = dvz_dy - dvy_dz
    curl_y = dvx_dz - dvz_dx
    curl_z = dvy_dx - dvx_dy

    return _apply_vector_mode(curl_x, curl_y, curl_z, grid_x, grid_y, grid_z, box_center, mode)


def compute_frame_curl(
    oxygen_positions,
    dipole_vectors,
    grid_x,
    grid_y,
    grid_z,
    dx,
    dy,
    dz,
    box_center,
    mode="magnitude",
):
    """Backward-compatible wrapper to compute the dipole curl for a single frame."""
    return compute_frame_property(
        oxygen_positions=oxygen_positions,
        dipole_vectors=dipole_vectors,
        grid_x=grid_x,
        grid_y=grid_y,
        grid_z=grid_z,
        dx=dx,
        dy=dy,
        dz=dz,
        box_center=box_center,
        property_name="curl",
        mode=mode,
    )

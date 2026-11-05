"""Multiprocessing worker for trajectory frame chunk processing."""

import time
import numpy as np

from dipole_calculus.core.grid import compute_frame_property
from dipole_calculus.io.exporters import export_cube


def process_chunk_worker(args):
    """Load an independent trajectory instance and calculate the chosen dipole property for a chunk of frames."""
    import MDAnalysis as mda
    (
        chunk_frames,
        reference_file,
        xtc_file,
        charges_array,
        center_sel,
        target_sel,
        grid_x,
        grid_y,
        grid_z,
        dx,
        dy,
        dz,
        min_coords,
        max_coords,
        box_center,
        grid_resolution,
        box_edges,
        property_name,
        mode,
        cube_interval,
    ) = args

    min_x, min_y, min_z = min_coords
    max_x, max_y, max_z = max_coords

    # Create an isolated MDAnalysis Universe instance for this worker process
    local_u = mda.Universe(reference_file, xtc_file)
    local_u.add_TopologyAttr("charges", charges_array)
    selection = local_u.select_atoms(target_sel)
    n_selection = selection.n_residues
    water_charges_matrix = selection.charges.reshape(n_selection, 3)

    results = []

    # Initialize 4D accumulator for vector mode on vector operators; 3D for scalar/magnitude modes
    is_vector_accumulation = (mode == "vector") and (property_name != "divergence")
    if is_vector_accumulation:
        chunk_prop_sum = np.zeros((3, grid_resolution[0], grid_resolution[1], grid_resolution[2]))
    else:
        chunk_prop_sum = np.zeros((grid_resolution[0], grid_resolution[1], grid_resolution[2]))

    for frame_idx in chunk_frames:
        frame_start_time = time.time()
        local_u.trajectory[frame_idx]

        # Compute molecular electric dipoles and oxygen coordinates for the frame
        selection_positions = selection.positions.reshape(n_selection, 3, 3)
        selection_dipoles = np.sum(selection_positions * water_charges_matrix[:, :, None], axis=1)
        oxygen_positions = selection_positions[:, 0, :]

        # Count molecules inside the spatial grid boundaries
        mask_x = (oxygen_positions[:, 0] >= min_x) & (oxygen_positions[:, 0] <= max_x)
        mask_y = (oxygen_positions[:, 1] >= min_y) & (oxygen_positions[:, 1] <= max_y)
        mask_z = (oxygen_positions[:, 2] >= min_z) & (oxygen_positions[:, 2] <= max_z)
        molecules_in_grid = int(np.sum(mask_x & mask_y & mask_z))

        # Interpolate dipoles onto grid and compute selected differential property
        frame_mean, frame_contribution, export_grid = compute_frame_property(
            oxygen_positions=oxygen_positions,
            dipole_vectors=selection_dipoles,
            grid_x=grid_x,
            grid_y=grid_y,
            grid_z=grid_z,
            dx=dx,
            dy=dy,
            dz=dz,
            box_center=box_center,
            property_name=property_name,
            mode=mode,
        )

        chunk_prop_sum += frame_contribution

        # Export per-frame CUBE file if requested at this interval
        if cube_interval and ((frame_idx + 1) % cube_interval == 0):
            cube_filename = f"frame_{frame_idx + 1}_{property_name}.cube"
            export_cube(cube_filename, export_grid, min_coords, max_coords, grid_resolution, box_edges)

        frame_elapsed = time.time() - frame_start_time
        results.append((frame_idx, frame_mean, molecules_in_grid, frame_elapsed))

        print(
            f"Frame: {frame_idx + 1:<6} | Waters: {molecules_in_grid:<6} | Time: {frame_elapsed:.4f}s",
            flush=True,
        )

    return results, chunk_prop_sum

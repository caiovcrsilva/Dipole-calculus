"""Main orchestrator for parallel electric dipole differential calculations across trajectories."""

from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import numpy as np

from dipole_calculus.core.grid import build_spatial_grid
from dipole_calculus.core.worker import process_chunk_worker
from dipole_calculus.io.exporters import export_cube, export_tcl_script
from dipole_calculus.io.topology import build_charges_from_top


def calculate_grid_property(
    reference_file,
    xtc_file,
    top_file,
    num_workers,
    box_edges,
    center_sel,
    target_sel,
    grid_resolution=(15, 15, 15),
    property_name="curl",
    mode="magnitude",
    cube_per_frame=None,
    generate_tcl=False,
):
    """Calculate a differential spatial property of water electric dipoles over time using parallel workers."""
    import MDAnalysis as mda
    print(f"Loading reference from: {reference_file}...", flush=True)
    load_start_time = time.time()
    # Load main universe to determine reference center of mass and total frames
    universe = mda.Universe(reference_file, xtc_file)
    load_end_time = time.time()
    print(f"-> Loaded in {load_end_time - load_start_time:.2f} seconds.", flush=True)

    charges = build_charges_from_top(top_file)
    total_frames = len(universe.trajectory)

    reference_atoms = universe.select_atoms(center_sel)
    if len(reference_atoms) > 0:
        box_center = reference_atoms.center_of_mass()
    else:
        print(f"WARNING: Selection '{center_sel}' is empty. Using system center of mass instead.", flush=True)
        box_center = universe.atoms.center_of_mass()

    # Generate 3D meshgrid and bounding box coordinates
    grid_x, grid_y, grid_z, dx, dy, dz, min_coords, max_coords = build_spatial_grid(
        box_center=box_center,
        box_edges=box_edges,
        grid_resolution=grid_resolution,
    )

    is_vector_accumulation = (mode == "vector") and (property_name != "divergence")
    effective_mode_label = "scalar" if property_name == "divergence" else mode

    print(
        f"\nRunning {property_name} ({effective_mode_label}) calculation for {total_frames} frames...",
        flush=True,
    )
    print(f"Grid: {grid_resolution[0]}x{grid_resolution[1]}x{grid_resolution[2]} points.", flush=True)

    # Split trajectory frames into chunks across available CPU workers
    all_frames = list(range(total_frames))
    chunk_size = max(1, total_frames // num_workers)
    chunks = [all_frames[i : i + chunk_size] for i in range(0, total_frames, chunk_size)]

    tasks = [
        (
            chunk,
            reference_file,
            xtc_file,
            charges,
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
            cube_per_frame,
        )
        for chunk in chunks
    ]

    print(f"\nInitializing parallel calculation using {num_workers} threads...", flush=True)
    traj_start_time = time.time()

    if is_vector_accumulation:
        total_prop_sum = np.zeros((3, grid_resolution[0], grid_resolution[1], grid_resolution[2]))
    else:
        total_prop_sum = np.zeros((grid_resolution[0], grid_resolution[1], grid_resolution[2]))

    unordered_results = []
    completed_frames = 0

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_chunk_worker, task) for task in tasks]

        for future in as_completed(futures):
            chunk_results, chunk_sum = future.result()
            unordered_results.extend(chunk_results)
            total_prop_sum += chunk_sum
            completed_frames += len(chunk_results)

    traj_total_time = time.time() - traj_start_time
    print(f"\nCalculation done! Total time: {traj_total_time:.2f} seconds.", flush=True)

    sorted_results = sorted(unordered_results, key=lambda x: x[0])

    # Compute spatio-temporal mean grid according to the chosen property and mode
    if is_vector_accumulation:
        mean_vec_grid = total_prop_sum / completed_frames
        final_export_grid = np.sqrt(mean_vec_grid[0] ** 2 + mean_vec_grid[1] ** 2 + mean_vec_grid[2] ** 2)
    else:
        final_export_grid = total_prop_sum / completed_frames

    mean_cube_filename = f"mean_{property_name}.cube"
    print(f"\nGenerating average CUBE file ({mean_cube_filename})...", flush=True)
    cube_start_time = time.time()
    export_cube(mean_cube_filename, final_export_grid, min_coords, max_coords, grid_resolution, box_edges)
    cube_end_time = time.time()
    print(f"CUBE file {mean_cube_filename} generated in {round(cube_end_time - cube_start_time, 2)}s", flush=True)

    if generate_tcl:
        tcl_filename = f"plot_{property_name}.vmd"
        print(f"Generating VMD TCL script ({tcl_filename})...", flush=True)
        export_tcl_script(mean_cube_filename, final_export_grid, tcl_filename=tcl_filename)

    return sorted_results


def calculate_grid_curl(
    reference_file,
    xtc_file,
    top_file,
    num_workers,
    box_edges,
    center_sel,
    target_sel,
    grid_resolution=(15, 15, 15),
    mode="magnitude",
    cube_per_frame=None,
    generate_tcl=False,
    property_name="curl",
):
    """Backward-compatible wrapper around calculate_grid_property."""
    return calculate_grid_property(
        reference_file=reference_file,
        xtc_file=xtc_file,
        top_file=top_file,
        num_workers=num_workers,
        box_edges=box_edges,
        center_sel=center_sel,
        target_sel=target_sel,
        grid_resolution=grid_resolution,
        property_name=property_name,
        mode=mode,
        cube_per_frame=cube_per_frame,
        generate_tcl=generate_tcl,
    )

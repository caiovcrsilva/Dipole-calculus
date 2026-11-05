"""Command-line interface and output formatting for dipole_calculus."""

import argparse
import os
import sys

from dipole_calculus.core.calculator import calculate_grid_property


def build_parser():
    """Create and configure the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Calculate spatial differential properties from water electric dipoles in a GROMACS trajectory."
    )

    parser.add_argument("-top", "--topol", required=True, help="Path for the global topology file (.top)")
    parser.add_argument("-ref", "--reference", required=True, help="Path for the reference coordinate file (.gro)")
    parser.add_argument("-traj", "--trajectory", required=True, help="Path for the trajectory file (.xtc)")
    parser.add_argument(
        "-g", "--grid", type=int, nargs=3, default=[80, 80, 80], help="Grid points along X Y Z"
    )
    parser.add_argument(
        "-cpus", "--parallel", type=int, default=os.cpu_count(), help="CPUs/Threads amount for parallel runs"
    )
    parser.add_argument(
        "-box", "--box", type=float, nargs=3, required=True, help="Total box size along X Y Z in Angstroms"
    )
    parser.add_argument(
        "-box_ref", "--box_ref", type=str, default="protein", help="Selection string for the grid center"
    )
    parser.add_argument(
        "-sel", "--sel", type=str, default="resname SOL", help="Selection string for the target molecules"
    )
    parser.add_argument("--gen_tcl", action="store_true", help="Generate VMD TCL script for visualization")
    parser.add_argument(
        "--property",
        type=str,
        choices=["curl", "divergence", "gradient", "laplacian"],
        default="curl",
        help="Differential property to compute (curl, divergence, gradient, or laplacian)",
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["magnitude", "vector", "dot_product"],
        default="magnitude",
        help="Calculation mode for vector properties (magnitude, vector, or dot_product)",
    )
    parser.add_argument(
        "--cube_per_frame", type=int, default=None, help="Frame interval to export individual CUBE files"
    )

    return parser


def print_results(sorted_results, mode, property_name="curl"):
    """Print the formatted summary table of per-frame results according to property and mode."""
    prop_title = property_name.capitalize()

    if property_name == "divergence":
        print(f"\n=== Final results ({prop_title} - Scalar Field) ===", flush=True)
        print(f"{'Frame':<10} | {'Average Divergence':<25} | {'Waters':<15} ", flush=True)
        print("-" * 55, flush=True)
        for r in sorted_results:
            print(f"{r[0] + 1:<10} | {r[1]:<25.6f} | {r[2]:<15} ", flush=True)
        return

    if mode == "vector":
        col_title = f"Average {prop_title} (X, Y, Z)"
        print(f"\n=== Final results ({prop_title} - Vector Mode) ===", flush=True)
        print(f"{'Frame':<10} | {col_title:<45} | {'Waters':<15} ", flush=True)
        print("-" * 75, flush=True)
        for r in sorted_results:
            frame_num = r[0] + 1
            x, y, z = r[1]
            vec_str = f"X: {x:8.5f}  Y: {y:8.5f}  Z: {z:8.5f}"
            print(f"{frame_num:<10} | {vec_str:<45} | {r[2]:<15} ", flush=True)
    elif mode == "dot_product":
        col_title = f"Average {prop_title} Dot Product"
        print(f"\n=== Final results ({prop_title} - Dot Product Mode) ===", flush=True)
        print(f"{'Frame':<10} | {col_title:<28} | {'Waters':<15} ", flush=True)
        print("-" * 58, flush=True)
        for r in sorted_results:
            print(f"{r[0] + 1:<10} | {r[1]:<28.6f} | {r[2]:<15} ", flush=True)
    else:
        col_title = f"Average {prop_title} Magnitude"
        print(f"\n=== Final results ({prop_title} - Magnitude Mode) ===", flush=True)
        print(f"{'Frame':<10} | {col_title:<28} | {'Waters':<15} ", flush=True)
        print("-" * 58, flush=True)
        for r in sorted_results:
            print(f"{r[0] + 1:<10} | {r[1]:<28.6f} | {r[2]:<15} ", flush=True)


def main(argv=None):
    """Parse CLI arguments and execute the dipole differential calculation pipeline."""
    raw_argv = sys.argv[1:] if argv is None else list(argv)
    parser = build_parser()
    args = parser.parse_args(raw_argv)

    # Warn without stopping if --mode is used with the scalar divergence property
    mode_explicitly_passed = any(arg == "--mode" or arg.startswith("--mode=") for arg in raw_argv)
    if args.property == "divergence" and (mode_explicitly_passed or args.mode in ("vector", "dot_product")):
        print(
            f"WARNING: Divergence is inherently a scalar field. The '--mode {args.mode}' "
            "option does not apply and will be ignored. Proceeding with standard scalar calculation.",
            flush=True,
        )

    try:
        sorted_results = calculate_grid_property(
            reference_file=args.reference,
            xtc_file=args.trajectory,
            top_file=args.topol,
            num_workers=args.parallel,
            box_edges=args.box,
            center_sel=args.box_ref,
            target_sel=args.sel,
            grid_resolution=args.grid,
            property_name=args.property,
            mode=args.mode,
            cube_per_frame=args.cube_per_frame,
            generate_tcl=args.gen_tcl,
        )
        print_results(sorted_results, args.mode, property_name=args.property)
    except Exception as exc:
        print(f"Error during execution: {exc}", flush=True)


if __name__ == "__main__":
    main()

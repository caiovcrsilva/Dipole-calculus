"""Dipole Calculus: Electric dipole differential operator calculations for GROMACS trajectories."""

from dipole_calculus.core.calculator import calculate_grid_curl, calculate_grid_property
from dipole_calculus.io.exporters import export_cube, export_tcl_script
from dipole_calculus.io.topology import build_charges_from_top, extract_charges_from_file

__version__ = "0.2.0"

__all__ = [
    "calculate_grid_curl",
    "calculate_grid_property",
    "build_charges_from_top",
    "extract_charges_from_file",
    "export_cube",
    "export_tcl_script",
]

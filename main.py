"""Entry-point script to run dipole_calculus directly from the repository root."""

import os
import sys

# Ensure the local 'src' directory is on the Python path when running without pip install
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from dipole_calculus.cli import main

if __name__ == "__main__":
    main()

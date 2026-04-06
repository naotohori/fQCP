"""
fQCP — Fast QCP RMSD and superposition.

Automatically uses the F2PY/Fortran backend if compiled (.so files present),
otherwise falls back to the pure Python/NumPy implementation.

Usage:
    from fQCP import calc_rmsd, calc_rotation, superimpose
"""

try:
    from ._fortran import calc_rmsd, calc_rotation, superimpose
    BACKEND = 'fortran'
except ImportError:
    from .qcp import calc_rmsd, calc_rotation, superimpose
    BACKEND = 'numpy'

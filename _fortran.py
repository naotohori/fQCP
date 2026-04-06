"""Thin wrappers around f2py modules, adapting to (N, 3) interface."""

import numpy as np
from .CalcRMSD import calcrmsd as _calcrmsd
from .CalcROT import calcrotation as _calcrotation
from .Superimpose import superimpose as _superimpose


def _to_fortran(coords):
    """Convert (N,3) array to (3,N) Fortran-contiguous copy."""
    return np.array(np.asarray(coords, dtype=np.float64).T, order='F')


def calc_rmsd(coords1, coords2):
    return _calcrmsd(_to_fortran(coords1), _to_fortran(coords2))


def calc_rotation(coords1, coords2):
    return _calcrotation(_to_fortran(coords1), _to_fortran(coords2))


def superimpose(coords1, coords2):
    c2 = _to_fortran(coords2)
    rmsd = _superimpose(_to_fortran(coords1), c2)
    return rmsd, c2.T.copy()

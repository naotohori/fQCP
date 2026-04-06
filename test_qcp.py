#!/usr/bin/env python
"""
Validation test for the QCP implementation.
Uses the same 7-atom test data as main.F90.
Expected RMSD: 0.719106 A

Run from the parent directory:
    python -m fQCP.test_qcp
Or directly:
    cd fQCP && python test_qcp.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import fQCP
from fQCP import calc_rmsd, calc_rotation, superimpose

print("Backend:", fQCP.BACKEND)
print()

# Test data from main.F90
frag_a = np.array([
    [-2.803, -15.373,  24.556],
    [ 0.893, -16.062,  25.147],
    [ 1.368, -12.371,  25.885],
    [-1.651, -12.153,  28.177],
    [-0.440, -15.218,  30.068],
    [ 2.551, -13.273,  31.372],
    [ 0.105, -11.330,  33.567],
])

frag_b = np.array([
    [-14.739, -18.673,  15.040],
    [-12.473, -15.810,  16.074],
    [-14.802, -13.307,  14.408],
    [-17.782, -14.852,  16.171],
    [-16.124, -14.617,  19.584],
    [-15.029, -11.037,  18.902],
    [-18.577, -10.001,  17.996],
])

EXPECTED_RMSD = 0.719106
TOL = 1e-4

print("RMSD should be {:.6f} A".format(EXPECTED_RMSD))
print()

# --- Test calc_rmsd ---
print("### Test calc_rmsd ###")
rmsd = calc_rmsd(frag_a, frag_b)
print("QCP rmsd: {:.6f}".format(rmsd))
assert abs(rmsd - EXPECTED_RMSD) < TOL, "RMSD mismatch: {:.6f}".format(rmsd)
print("PASS")
print()

# --- Test calc_rotation ---
print("### Test calc_rotation ###")
rmsd, mat = calc_rotation(frag_a, frag_b)
print("QCP rmsd: {:.6f}".format(rmsd))
assert abs(rmsd - EXPECTED_RMSD) < TOL, "RMSD mismatch: {:.6f}".format(rmsd)
print("Rotation-translation matrix (4x4):")
for i in range(4):
    print("  " + "  ".join("{:12.8f}".format(mat[i, j]) for j in range(4)))

# Verify by applying the matrix
transformed = frag_b @ mat[:3, :3].T + mat[:3, 3]
euc_dist = np.sqrt(np.sum((frag_a - transformed)**2) / len(frag_a))
print("Explicit RMSD from transformed coords: {:.6f}".format(euc_dist))
assert abs(euc_dist - EXPECTED_RMSD) < TOL
print("PASS")
print()

# --- Test superimpose ---
print("### Test superimpose ###")
rmsd, frag_b_tr = superimpose(frag_a, frag_b)
print("QCP rmsd: {:.6f}".format(rmsd))
assert abs(rmsd - EXPECTED_RMSD) < TOL, "RMSD mismatch: {:.6f}".format(rmsd)

print("Rotation matrix (3x3 block from calc_rotation):")
for i in range(3):
    print("  " + "  ".join("{:12.8f}".format(mat[i, j]) for j in range(3)))

print("Coords 2 after superimpose:")
for i in range(len(frag_b_tr)):
    print("  {:8.3f} {:8.3f} {:8.3f}".format(*frag_b_tr[i]))

euc_dist = np.sqrt(np.sum((frag_a - frag_b_tr)**2) / len(frag_a))
print("Explicit RMSD from superimposed coords: {:.6f}".format(euc_dist))
assert abs(euc_dist - EXPECTED_RMSD) < TOL
print("PASS")
print()

print("All tests passed.")

#!/usr/bin/env python
"""Benchmark: pure Python/NumPy qcp vs f2py Fortran backend."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import time
import numpy as np

from fQCP.qcp import calc_rotation as calc_rotation_numpy

try:
    from fQCP._fortran import calc_rotation as calc_rotation_fortran
    HAS_F2PY = True
except ImportError:
    HAS_F2PY = False
    print("WARNING: f2py modules not found, skipping Fortran benchmark.")
    print("         Run 'make f2py' to build them.\n")


def make_random_coords(n_atoms, seed=42):
    rng = np.random.default_rng(seed)
    c1 = rng.standard_normal((n_atoms, 3)) * 10.0
    angle = 0.3
    R = np.array([[np.cos(angle), -np.sin(angle), 0],
                  [np.sin(angle),  np.cos(angle), 0],
                  [0, 0, 1]])
    c2 = c1 @ R.T + rng.standard_normal(3) * 2.0 + rng.standard_normal((n_atoms, 3)) * 0.5
    return c1, c2


def bench(func, args, n_repeat):
    func(*args)  # warmup
    t0 = time.perf_counter()
    for _ in range(n_repeat):
        func(*args)
    return (time.perf_counter() - t0) / n_repeat


print("{:>8s}  {:>6s}  {:>12s}  {:>12s}  {:>8s}".format(
    "N_atoms", "calls", "NumPy (ms)", "Fortran (ms)", "ratio"))
print("-" * 58)

for n_atoms in [7, 50, 200, 1000, 5000]:
    c1, c2 = make_random_coords(n_atoms)
    n_repeat = max(100, 10000 // n_atoms)

    t_py = bench(calc_rotation_numpy, (c1, c2), n_repeat) * 1000

    if HAS_F2PY:
        t_f = bench(calc_rotation_fortran, (c1, c2), n_repeat) * 1000
        ratio = t_py / t_f
        print("{:8d}  {:6d}  {:12.4f}  {:12.4f}  {:7.1f}x".format(
            n_atoms, n_repeat, t_py, t_f, ratio))
    else:
        print("{:8d}  {:6d}  {:12.4f}  {:>12s}  {:>8s}".format(
            n_atoms, n_repeat, t_py, "N/A", "N/A"))

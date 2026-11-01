#!/usr/bin/env python
"""Compare the mean/rms maps saved by two runs of bench_rms_mean_map.py

Prints, for every map, whether the non-finite entries agree, whether the maps
are bit-identical, the largest difference relative to the local rms, the
largest relative difference, and the boxes whose rms differs by more than
--rtol (default 1e-4).

Usage: python compare_maps.py <reference.npz> <new.npz> [--rtol 1e-4]
"""
import sys

import numpy as np

ref = np.load(sys.argv[1])
new = np.load(sys.argv[2])
rtol = float(sys.argv[sys.argv.index('--rtol') + 1]) if '--rtol' in sys.argv else 1e-4

for key in ref.files:
    x = ref[key].astype(float)
    y = new[key].astype(float)
    rms = ref[key.replace('_mean', '_rms')].astype(float)
    fin = np.isfinite(x)
    diff = np.abs(x - y)[fin]
    print(f'{key:18s} shape={x.shape} '
          f'nonfinite_same={np.array_equal(np.isfinite(x), np.isfinite(y))} '
          f'identical={np.array_equal(ref[key], new[key], equal_nan=True)} '
          f'max|diff|/rms={np.max(diff / np.abs(rms[fin])):.1e} '
          f'max_rel={np.max(diff / np.maximum(np.abs(x[fin]), 1e-30)):.1e}')
    if key.endswith('_rms'):
        rel = np.abs(x - y) / np.abs(x)
        for idx in np.argwhere(rel > rtol):
            idx = tuple(idx)
            print(f'    box {idx}: {x[idx]:.6g} -> {y[idx]:.6g} (rel {rel[idx]:.1e})')

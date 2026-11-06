#!/usr/bin/env python
"""Replay recorded lmder_fit calls with the C++ `_cbdsm` and with
`bdsf/mgfunction.py`, and compare their results and run times.

Every recorded call is replayed with a new MGFunction, so the cache of fit
results of the Python version is not used; this compares the speed of single
fits.

Usage: replay_fits.py <records.pkl> <path to mgfunction.py> [repeats]

Run with a Python environment that has a PyBDSF version with `_cbdsm`.
"""
import gc
import importlib.util
import pickle
import sys
import time

import numpy as np

import bdsf._cbdsm as C

spec = importlib.util.spec_from_file_location("mgfunction", sys.argv[2])
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)
with open(sys.argv[1], 'rb') as f:
    records = pickle.load(f)
repeats = int(sys.argv[3]) if len(sys.argv) > 3 else 3


def run(M, r):
    fcn = M.MGFunction(r['data'], r['mask'], 1)
    fcn.find_peak()
    for gtype, p in zip(r['gtypes'], r['p0']):
        fcn.add_gaussian(M.Gtype(gtype), p)
    ok = M.lmder_fit(fcn, final=r['final'], verbose=0)
    return ok, fcn.parameters, fcn.find_peak()


def best_times(M):
    best = np.full(len(records), np.inf)
    for _ in range(repeats):
        for i, r in enumerate(records):
            start = time.perf_counter()
            run(M, r)
            best[i] = min(best[i], time.perf_counter() - start)
    return best


gc.disable()
identical = 0
status_differs = 0
max_rel_diff = []
for r in records:
    ok_c, par_c, peak_c = run(C, r)
    ok_p, par_p, peak_p = run(P, r)
    identical += (ok_c == ok_p and par_c == par_p and peak_c == peak_p)
    status_differs += (ok_c != ok_p)
    par_c, par_p = np.array(par_c), np.array(par_p)
    max_rel_diff.append(np.max(np.abs(par_c - par_p) / np.maximum(np.abs(par_c), 1e-12))
                        if par_c.size else 0.0)
max_rel_diff = np.array(max_rel_diff)
print(f"{len(records)} fits; bit-identical results: {identical}; "
      f"different convergence status: {status_differs}; "
      f"max rel param diff > 1e-10: {(max_rel_diff > 1e-10).sum()}")

t_c = best_times(C)
t_p = best_times(P)
print(f"total time  C++: {t_c.sum():.3f} s   Python: {t_p.sum():.3f} s   "
      f"ratio Python/C++: {t_p.sum() / t_c.sum():.2f}")
npix = np.array([(~r['mask']).sum() for r in records])
for lo, hi in ((0, 50), (50, 200), (200, 1000), (1000, np.inf)):
    sel = (npix >= lo) & (npix < hi)
    if sel.any():
        print(f"  pixels [{lo}, {hi}): {sel.sum():6d} fits  C++ {t_c[sel].sum():7.3f} s  "
              f"Python {t_p[sel].sum():7.3f} s  ratio {t_p[sel].sum() / t_c[sel].sum():.2f}")

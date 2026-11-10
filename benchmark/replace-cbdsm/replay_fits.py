#!/usr/bin/env python
"""Replay recorded lmder_fit calls with two implementations of MGFunction and
lmder_fit, and compare their results, run times and memory use.

Every recorded call is replayed with a new MGFunction, so the cache of fit
results of the Python version is not used; this compares single fits.

Usage: replay_fits.py <records.pkl> <impl A> <impl B> [repeats]

Each implementation is the path of a Python file (e.g. bdsf/mgfunction.py, or
an older version from `git show`), or the name of a module (e.g. bdsf._cbdsm,
which needs a PyBDSF version that still has it).
"""
import gc
import importlib
import importlib.util
import pickle
import sys
import time
import tracemalloc

import numpy as np


def load(impl, name):
    if not impl.endswith('.py'):
        return importlib.import_module(impl)
    spec = importlib.util.spec_from_file_location(name, impl)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


with open(sys.argv[1], 'rb') as f:
    records = pickle.load(f)
A, B = load(sys.argv[2], 'impl_a'), load(sys.argv[3], 'impl_b')
repeats = int(sys.argv[4]) if len(sys.argv) > 4 else 3


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


identical = status_differs = 0
max_rel_diff = []
for r in records:
    ok_a, par_a, peak_a = run(A, r)
    ok_b, par_b, peak_b = run(B, r)
    # repr() is exact for floats, and also equal for NaNs
    identical += repr((ok_a, par_a, peak_a)) == repr((ok_b, par_b, peak_b))
    status_differs += ok_a != ok_b
    par_a, par_b = np.array(par_a), np.array(par_b)
    with np.errstate(invalid='ignore'):
        diff = np.abs(par_a - par_b) / np.maximum(np.abs(par_a), 1e-12)
    diff[(par_a == par_b) | (np.isnan(par_a) & np.isnan(par_b))] = 0.0
    max_rel_diff.append(np.nan_to_num(diff, nan=np.inf).max() if diff.size else 0.0)
print(f"{len(records)} fits; bit-identical results: {identical}; "
      f"different convergence status: {status_differs}; "
      f"max rel param diff > 1e-10: {(np.array(max_rel_diff) > 1e-10).sum()}")

gc.disable()
t_a, t_b = best_times(A), best_times(B)
gc.enable()
print(f"total time  A: {t_a.sum():.3f} s   B: {t_b.sum():.3f} s   "
      f"ratio B/A: {t_b.sum() / t_a.sum():.3f}")
npix = np.array([(~r['mask']).sum() for r in records])
for lo, hi in ((0, 50), (50, 200), (200, 1000), (1000, np.inf)):
    sel = (npix >= lo) & (npix < hi)
    if sel.any():
        print(f"  pixels [{lo}, {hi}): {sel.sum():6d} fits  A {t_a[sel].sum():7.3f} s  "
              f"B {t_b[sel].sum():7.3f} s  ratio {t_b[sel].sum() / t_a[sel].sum():.3f}")

# Peak memory (traced by Python, so not that of C++) of the largest fits
largest = np.argsort(-npix * [len(r['gtypes']) for r in records])[:5]
for M, name in ((A, 'A'), (B, 'B')):
    peaks = []
    for i in largest:
        tracemalloc.start()
        run(M, records[i])
        peaks.append(tracemalloc.get_traced_memory()[1] // 1024)
        tracemalloc.stop()
    print(f"peak traced memory of the 5 largest fits, {name}: {peaks} KiB")

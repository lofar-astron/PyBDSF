#!/usr/bin/env python
"""Compare speed and results of two versions of functions.bstat

Runs bstat from both packages on a sample of 150x150 boxes (every 25th row of
boxes, step 15) and on the full image with a mask, and reports the time per
box, whether the iteration counts agree, and the relative differences of
m_raw, r_raw, m (clipped mean) and r (clipped rms).

Usage:
  python -I compare_bstat.py <pkg_orig dir> <pkg_new dir> <fits file>
"""
import importlib.util
import os
import sys
import time

import numpy as np
from astropy.io import fits


def load_functions(pkg_dir):
    """Load bdsf/functions.py from pkg_dir as a stand-alone module"""
    path = os.path.join(pkg_dir, 'bdsf', 'functions.py')
    spec = importlib.util.spec_from_file_location('functions_' + str(abs(hash(pkg_dir))), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


f_orig = load_functions(sys.argv[1])
f_new = load_functions(sys.argv[2])
data = fits.getdata(sys.argv[3])[0, 0].astype(np.float32)
kappa = -30.0

inds = [(i*15, i*15+150, j*15, j*15+150)
        for i in range(0, (data.shape[0]-150)//15, 25)
        for j in range((data.shape[1]-150)//15)]
results = {}
for name, f in [('orig', f_orig.bstat), ('new', f_new.bstat)]:
    t = time.time()
    results[name] = np.array([f(data[a:b, c:d], None, kappa) for a, b, c, d in inds])
    print(f'{name}: {(time.time() - t) / len(inds) * 1e3:.3f} ms per box')

orig, new = results['orig'], results['new']
rel = np.abs(orig - new)[:, :4] / np.abs(orig[:, :4])
print('iteration counts equal:', np.array_equal(orig[:, 4], new[:, 4]))
print('max rel diff    (m_raw, r_raw, m, r):', rel.max(0))
print('median rel diff (m_raw, r_raw, m, r):', np.median(rel, 0))

mask = np.zeros(data.shape, bool)
mask[:3000, :2000] = True
for name, f in [('orig', f_orig.bstat), ('new', f_new.bstat)]:
    t = time.time()
    res = f(data, mask, kappa)
    print(f'{name} full image with mask: {time.time() - t:.2f} s', res)

#!/usr/bin/env python
"""Time Op_rmsimage.rms_mean_map() on a FITS image and save the maps

Cases:
  full_nomask : whole image, no mask, rms_box (150, 15) (extrapolation path)
  full_mask   : whole image with a masked (NaN) region, rms_box (60, 20)
  cut_pad     : 1200x1000 cutout, rms_box (150, 50) (padded-array path)
  tiled       : (with --tiled) image tiled 2x2, no mask, rms_box (150, 15)

Usage:
  python -I bench_rms_mean_map.py <pkg dir> <label> <fits file> <out.npz> [--tiled]

<pkg dir> is a directory containing a bdsf package (see setup_overlays.sh).
Compare two output files with compare_maps.py.
"""
import sys
import time

sys.path.insert(0, sys.argv[1])

import numpy as np
from astropy.io import fits
import bdsf
from bdsf.rmsimage import Op_rmsimage

label, fitsname, outname = sys.argv[2:5]
tiled = '--tiled' in sys.argv[5:]
print('Using', bdsf.__file__, file=sys.stderr)

data = fits.getdata(fitsname)[0, 0].astype(np.float32)
op = Op_rmsimage()
kappa = -30.0

if tiled:
    cases = {'tiled': (np.tile(data, (2, 2)), None, (150, 15))}
else:
    mask = np.zeros(data.shape, bool)
    mask[:400, :3000] = True
    mask[2000:2003, :] = True
    data_masked = data.copy()
    data_masked[mask] = np.nan
    cases = {
        'full_nomask': (data, None, (150, 15)),
        'full_mask': (data_masked, mask, (60, 20)),
        'cut_pad': (data[:1200, :1000].copy(), None, (150, 50)),
    }

out = {}
for name, (arr, msk, box) in cases.items():
    t = time.time()
    axes, mean_map, rms_map = op.rms_mean_map(arr, msk, kappa, box, None)
    print(f'{label} {name}: {time.time() - t:.1f} s', flush=True)
    out[name + '_mean'] = mean_map
    out[name + '_rms'] = rms_map
np.savez(outname, **out)

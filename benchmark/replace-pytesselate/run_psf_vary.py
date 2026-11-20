"""Run PyBDSF on NGC6048 with PSF variation, and report the time of each step.

Usage: run_psf_vary.py <image> [psf_itess_method]

The options are those of `benchmark/run_pybdsf.py`, with `psf_vary_do=True`.
`psf_itess_method` (default 0) selects the tessellation: 0 (simple, unity
weights), 1 (roundness), 2 (log10 weights) or 3 (sqrt(log10) weights).
"""
import sys
import timeit

import bdsf
from bdsf import psf_vary

image = sys.argv[1]
itess_method = int(sys.argv[2]) if len(sys.argv) > 2 else 0

# Time spent in the tessellation itself
tess_time = [0.0]
tesselate = psf_vary.Op_psf_vary.tesselate


def timed_tesselate(*args, **kwargs):
    t0 = timeit.default_timer()
    result = tesselate(*args, **kwargs)
    tess_time[0] += timeit.default_timer() - t0
    return result


psf_vary.Op_psf_vary.tesselate = timed_tesselate

start = timeit.default_timer()
img = bdsf.process_image(image,
                         thresh_isl=4.0,
                         thresh_pix=5.0,
                         rms_box=(150, 15),
                         rms_map=True,
                         mean_map='zero',
                         ini_method='intensity',
                         adaptive_rms_box=True,
                         adaptive_thresh=150,
                         rms_box_bright=(60, 15),
                         group_by_isl=False,
                         group_tol=10.0,
                         output_opts=True,
                         output_all=True,
                         atrous_do=True,
                         atrous_jmax=4,
                         flagging_opts=True,
                         flag_maxsize_fwhm=0.5,
                         advanced_opts=True,
                         blank_limit=None,
                         frequency=143.65e6,
                         psf_vary_do=True,
                         psf_itess_method=itess_method,
                         clobber=True)
total = timeit.default_timer() - start
print(f'PyBDSF running time: {total:.1f} s, of which tessellation: {tess_time[0]:.2f} s')

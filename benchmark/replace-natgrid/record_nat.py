"""Run benchmark/run_pybdsf.py, recording every natgrid call (inputs, output, CPU time).

Usage (with a PyBDSF installation that has `bdsf.nat`):
    python -I record_nat.py <output dir> <run dir containing run_pybdsf.py and the image>

The coordinates are stored as contiguous copies (what psf_vary meant to pass;
natgridmodule itself read the strided views wrongly, see DESCRIPTION.md).
"""
import os, sys, time, uuid, runpy
import numpy as N
from bdsf import nat
OUT = sys.argv[1]
orig_rgrd = nat.Natgrid.rgrd
def rgrd(self, z, *a, **k):
    xi, yi, xo, yo = self.xi.copy(), self.yi.copy(), self.xo.copy(), self.yo.copy()
    t0 = time.process_time()
    r = orig_rgrd(self, z, *a, **k)
    dt = time.process_time() - t0
    N.savez(os.path.join(OUT, uuid.uuid4().hex + '.npz'), xi=xi, yi=yi, xo=xo, yo=yo, z=z, out=r, t=dt)
    return r
nat.Natgrid.rgrd = rgrd
os.chdir(sys.argv[2])
sys.argv = ['run_pybdsf.py']
runpy.run_path('run_pybdsf.py', run_name='__main__')

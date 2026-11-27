"""Time and peak memory of one natgrid call (C extension or bdsf/natgrid.py of this tree).

Usage: python -I bench_call.py <c|py> <npz file from record_nat.py>
"""
import os, sys, time, resource, importlib.util
import numpy as N
which, f = sys.argv[1], sys.argv[2]
d = N.load(f)
xi, yi, z, xo, yo = (N.ascontiguousarray(d[k]) for k in ('xi', 'yi', 'z', 'xo', 'yo'))
if which == 'c':
    from bdsf import nat
    fn = lambda: nat.Natgrid(xi, yi, xo, yo).rgrd(z)
else:
    spec = importlib.util.spec_from_file_location('natgrid', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'bdsf', 'natgrid.py'))
    ng = importlib.util.module_from_spec(spec); spec.loader.exec_module(ng)
    fn = lambda: ng.natgrid(xi, yi, z, xo, yo)
time.sleep(2)  # let the OpenBLAS threads started by `import numpy` settle
base = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
t0 = time.process_time(); w0 = time.perf_counter()
out = fn()
t = time.process_time() - t0; w = time.perf_counter() - w0
peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print('%-2s cpu %.2f s  wall %.2f s  peak RSS increase %.0f MB (output %.0f MB)' % (which, t, w, (peak - base) / 1024, out.nbytes / 2**20))

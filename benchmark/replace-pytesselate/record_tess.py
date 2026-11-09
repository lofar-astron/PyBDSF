"""Record all calls to the Fortran extension `_pytesselate` during a PyBDSF run.

Usage: record_tess.py <out.pkl> <image> [process_image kwargs as a dict literal]

Each record is a dict with the function name, its arguments and its result
(the `volrank` array). The records are used by `bench_tess.py`, which needs
only the arguments, and checks the results itself.
"""
import ast
import pickle
import sys
import time

import bdsf
from bdsf import _pytesselate

out, image = sys.argv[1:3]
kwargs = ast.literal_eval(sys.argv[3]) if len(sys.argv) > 3 else {}
records = []


def recorder(name):
    func = getattr(_pytesselate, name)

    def wrapper(*args):
        t0 = time.perf_counter()
        result = func(*args)
        elapsed = time.perf_counter() - t0
        records.append(dict(func=name, args=args, result=result, time=elapsed))
        print(f'{name}: shape=({args[0]}, {args[1]}), ngens={len(args[2])}, '
              f'args={args[-2:]}, {elapsed:.3f} s', flush=True)
        return result
    return wrapper


# Both the Fortran names (`_pytess.pytess_*`) and those of bdsf.tessellate
# (`tessellate.tessellate_*`), so that this works with either version of
# psf_vary.py
bdsf.psf_vary._pytess = bdsf.psf_vary.tessellate = type('pytess', (), {
    prefix + kind: staticmethod(recorder('pytess_' + kind))
    for prefix in ('pytess_', 'tessellate_') for kind in ('simple', 'roundness')})
bdsf.process_image(image, quiet=True, psf_vary_do=True, **kwargs)
with open(out, 'wb') as f:
    pickle.dump(records, f)
print(f'Saved {len(records)} records to {out}')

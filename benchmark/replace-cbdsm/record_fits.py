#!/usr/bin/env python
"""Record all fit calls (lmder_fit or lm_fit_with_derivatives) of a PyBDSF run.

For each call, the data, mask, Gaussian types, initial parameters and
precision are saved, together with the result. This works with the C++
`_cbdsm` and with `bdsf.mgfunction`, whichever the installed PyBDSF uses.
The run uses a single core, so that all calls are recorded.

Usage: record_fits.py <out.pkl> <fits file> [process_image kwargs as a dict]

Run with `-I`, so that the installed `bdsf` is used, not the repository's.
"""
import ast
import pickle
import sys

import numpy as np

import bdsf

try:
    import bdsf._cbdsm as impl  # gausfit imports from it in each call
except ImportError:
    import bdsf.gausfit as impl  # gausfit imported from bdsf.mgfunction

# The fit function is called lmder_fit in _cbdsm and older versions of
# bdsf.mgfunction, and lm_fit_with_derivatives in newer versions
FIT_NAME = ('lmder_fit' if hasattr(impl, 'lmder_fit')
            else 'lm_fit_with_derivatives')
records = []
MGFunction, fit = impl.MGFunction, getattr(impl, FIT_NAME)


class RecordingMGFunction:
    """Wraps an MGFunction to keep track of its data and Gaussian types."""
    def __init__(self, data, mask, weight=1.0):
        self.fcn = MGFunction(data, mask, weight)
        self.data = np.array(data)
        self.mask = np.array(mask)
        self.gtypes = []

    def reset(self):
        self.gtypes = []
        self.fcn.reset()

    def add_gaussian(self, gtype, parameters):
        self.gtypes.append(int(gtype))
        self.fcn.add_gaussian(gtype, parameters)

    def __getattr__(self, name):
        return getattr(self.fcn, name)


def recording_fit(fcn, final=False, verbose=1):
    p0 = fcn.parameters
    ok = fit(fcn.fcn, final=final, verbose=verbose)
    records.append(dict(data=fcn.data, mask=fcn.mask, gtypes=list(fcn.gtypes),
                        p0=p0, final=bool(final), ok=bool(ok),
                        p1=fcn.parameters))
    return ok


impl.MGFunction = RecordingMGFunction
setattr(impl, FIT_NAME, recording_fit)

kwargs = ast.literal_eval(sys.argv[3]) if len(sys.argv) > 3 else {}
bdsf.process_image(sys.argv[2], ncores=1, quiet=True, **kwargs)
with open(sys.argv[1], 'wb') as f:
    pickle.dump(records, f, protocol=5)
print(len(records), FIT_NAME, 'calls recorded')

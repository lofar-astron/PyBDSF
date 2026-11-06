#!/usr/bin/env python
"""Record all lmder_fit calls of a PyBDSF run that uses the C++ `_cbdsm`.

For each call, the data, mask, Gaussian types, initial parameters and
precision are saved, together with the result. The run uses a single core, so
that all calls are recorded.

Usage: record_fits.py <out.pkl> <fits file> [process_image kwargs as a dict]

Run with a Python environment that has a PyBDSF version with `_cbdsm`, and
with `-I`, so that the repository's `bdsf` isn't imported.
"""
import ast
import pickle
import sys

import numpy as np

import bdsf
import bdsf._cbdsm as C

records = []
_MGFunction = C.MGFunction
_lmder_fit = C.lmder_fit


class RecordingMGFunction:
    """Wraps C.MGFunction to keep track of its data and Gaussian types."""
    def __init__(self, data, mask, weight=1.0):
        self.fcn = _MGFunction(data, mask, weight)
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

    def __len__(self):
        return len(self.fcn)


def recording_lmder_fit(fcn, final=False, verbose=1):
    p0 = fcn.fcn.parameters
    ok = _lmder_fit(fcn.fcn, final=final, verbose=verbose)
    records.append(dict(data=fcn.data, mask=fcn.mask, gtypes=list(fcn.gtypes),
                        p0=p0, final=bool(final), ok=bool(ok),
                        p1=fcn.fcn.parameters))
    return ok


C.MGFunction = RecordingMGFunction
C.lmder_fit = recording_lmder_fit

kwargs = ast.literal_eval(sys.argv[3]) if len(sys.argv) > 3 else {}
bdsf.process_image(sys.argv[2], ncores=1, quiet=True, **kwargs)
with open(sys.argv[1], 'wb') as f:
    pickle.dump(records, f, protocol=5)
print(len(records), 'lmder_fit calls recorded')

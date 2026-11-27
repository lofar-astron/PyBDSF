"""Create test/natgrid_reference.npz: outputs of the C extension natgridmodule.

Run with a PyBDSF installation that still has `bdsf.nat` (the C extension):

    python -I make_reference.py <output file>

Each case is stored as <name>_x, <name>_y, <name>_z, <name>_xo, <name>_yo
(inputs) and <name>_out (output of natgridd). All arrays are passed to the C
code as contiguous float64 arrays (natgridmodule does not check contiguity).
"""
import sys

import numpy as np
from bdsf import nat

rng = np.random.default_rng(20261124)
cases = {}


def add(name, x, y, z, xo, yo):
    x, y, z, xo, yo = (np.ascontiguousarray(a, dtype=np.float64) for a in (x, y, z, xo, yo))
    out = nat.Natgrid(x, y, xo, yo).rgrd(z)
    for key, val in zip(('x', 'y', 'z', 'xo', 'yo', 'out'), (x, y, z, xo, yo, out)):
        cases[name + '_' + key] = np.array(val)


# Random points inside the grid
x, y = rng.uniform(0, 100, (2, 20))
add('random', x, y, rng.normal(size=20), np.arange(0.0, 101.0, 4.0), np.arange(0.0, 101.0, 5.0))

# As in PyBDSF: few points (tile centers), grid much larger than their hull
x, y = rng.uniform(1500, 4500, (2, 13))
add('tiles', x, y, 2 + 0.3 * rng.normal(size=13), np.arange(0.0, 6000.0, 250.0), np.arange(0.0, 6000.0, 200.0))

# Points partly outside the grid
x, y = rng.uniform(-50, 150, (2, 30))
add('outside', x, y, np.sin(x / 30) + np.cos(y / 20), np.linspace(0, 100, 23), np.linspace(10, 90, 17))

# Very different extents in x and y (natgrid rescales the axes)
x = rng.uniform(0, 1000, 25)
y = rng.uniform(0, 60, 25)
add('aspect', x, y, x / 100 + np.sin(y / 10), np.arange(0.0, 1001.0, 40.0), np.arange(0.0, 61.0, 3.0))

# (Decreasing output coordinates are not included: natgridmodule crashes.)

# Grid points at the data points
xg, yg = np.arange(0.0, 50.0, 2.0), np.arange(0.0, 40.0, 2.0)
idx = rng.choice(xg.size * yg.size, 12, replace=False)
x, y = xg[idx // yg.size] + 0.0, yg[idx % yg.size] + 0.0
add('ondata', x, y, rng.normal(size=12), xg, yg)

np.savez_compressed(sys.argv[1], **cases)
print('Wrote %d cases to %s' % (len(cases) // 6, sys.argv[1]))

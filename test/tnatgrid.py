"""Tests for bdsf.natgrid, the replacement of the C extension natgridmodule.

The results are compared with outputs of the C extension, stored in
natgrid_reference.npz (created with benchmark/replace-natgrid/make_reference.py).
Natural neighbor interpolation also reproduces linear functions exactly, and
the interpolant is continuous; these properties are tested as well.

Run with: python tnatgrid.py (exit status 0 if all tests pass)
"""
import os
import sys

import numpy as np

from bdsf import natgrid as ng

HERE = os.path.dirname(os.path.abspath(__file__))
failures = []


def check(name, ok, detail=''):
    print('%-50s %s %s' % (name, 'ok' if ok else 'FAILED', detail))
    if not ok:
        failures.append(name)


def test_reference():
    """Compare with the outputs of the C extension."""
    ref = np.load(os.path.join(HERE, 'natgrid_reference.npz'))
    names = sorted({k.rsplit('_', 1)[0] for k in ref.files})
    for name in names:
        x, y, z, xo, yo, out = (ref[name + '_' + k] for k in ('x', 'y', 'z', 'xo', 'yo', 'out'))
        res = ng.natgrid(x, y, z, xo, yo)
        err = np.abs(res - out).max() / np.abs(out).max()
        check('reference: ' + name, res.shape == out.shape and err < 1e-9, '(max rel. diff %.1e)' % err)


def test_linear():
    """Linear functions are reproduced (also outside the convex hull)."""
    rng = np.random.default_rng(1)
    x, y = rng.uniform(100, 200, (2, 15))
    xo, yo = np.linspace(0, 300, 61), np.linspace(50, 250, 41)
    res = ng.natgrid(x, y, 3 - 0.02 * x + 0.05 * y, xo, yo)
    exp = 3 - 0.02 * xo[:, None] + 0.05 * yo[None, :]
    err = np.abs(res - exp).max()
    # The data coordinates are perturbed (as in natgrid) by up to about
    # 1e-5 * 300 / 2, which changes the values by up to about 1e-4
    check('linear function', err < 2e-4, '(max diff %.1e)' % err)


def test_multiple_values():
    """Several sets of values at once give the same as one by one."""
    rng = np.random.default_rng(2)
    x, y = rng.uniform(0, 50, (2, 10))
    z = rng.normal(size=(10, 3))
    xo, yo = np.arange(0.0, 50.0, 1.5), np.arange(0.0, 50.0, 2.5)
    res = ng.natgrid(x, y, z, xo, yo)
    ok = res.shape == (xo.size, yo.size, 3)
    for k in range(3):
        ok &= np.allclose(res[:, :, k], ng.natgrid(x, y, z[:, k], xo, yo), rtol=0, atol=1e-12)
    check('multiple sets of values', ok)


def test_degenerate():
    """Grid points on edges of the triangulation give finite, continuous values.

    Without the perturbation of the coordinates, data points on a lattice
    have edges through many grid points.
    """
    xd, yd = np.meshgrid(np.arange(0.0, 50.0, 10.0), np.arange(0.0, 50.0, 10.0), indexing='ij')
    pts = np.column_stack((xd.ravel(), yd.ravel()))
    vals = np.sin(pts[:, :1] / 15) + np.cos(pts[:, 1:] / 25)
    big = np.array([[-1000.0, -1000.0], [-1000.0, 5000.0], [5000.0, -1000.0]])
    pts = np.vstack((pts, big))
    vals = np.vstack((vals, np.zeros((3, 1))))
    xo = yo = np.arange(0.0, 41.0, 5.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        res = ng._interpolate(pts, vals, xo, yo)[:, :, 0]
        near = ng._interpolate(pts, vals, xo + 1e-6, yo + 0.4e-6)[:, :, 0]
    ok = np.isfinite(res).all() and np.abs(res - near).max() < 1e-4
    check('grid points on edges', ok, '(max diff with shifted grid %.1e)' % np.abs(res - near).max())


def test_noncontiguous():
    """Strided input arrays (as from numpy.transpose) are handled correctly."""
    rng = np.random.default_rng(4)
    coords = rng.uniform(0, 100, (12, 2))
    z = rng.normal(size=12)
    xo, yo = np.arange(0.0, 100.0, 7.0), np.arange(0.0, 100.0, 9.0)
    x, y = np.transpose(coords)
    res = ng.natgrid(x, y, z, xo, yo)
    exp = ng.natgrid(x.copy(), y.copy(), z, xo, yo)
    check('non-contiguous input', np.array_equal(res, exp))


def test_errors():
    """Invalid input raises ValueError."""
    xo = yo = np.arange(0.0, 10.0)
    for name, args in (('too few points', ([1.0, 2, 3], [1.0, 5, 2], [1.0, 2, 3], xo, yo)),
                       ('duplicate points', ([1.0, 2, 3, 2], [1.0, 5, 2, 5], [1.0, 2, 3, 4], xo, yo)),
                       ('non-monotonic grid', ([1.0, 2, 3, 7], [1.0, 5, 2, 8], [1.0, 2, 3, 4], [0.0, 2, 1], yo)),
                       ('mismatched lengths', ([1.0, 2, 3, 7], [1.0, 5, 2, 8], [1.0, 2, 3], xo, yo))):
        try:
            ng.natgrid(*args)
            ok = False
        except ValueError:
            ok = True
        check('error: ' + name, ok)


if __name__ == '__main__':
    test_reference()
    test_linear()
    test_multiple_values()
    test_degenerate()
    test_noncontiguous()
    test_errors()
    if failures:
        print('%d test(s) FAILED' % len(failures))
        sys.exit(1)
    print('All tests passed')

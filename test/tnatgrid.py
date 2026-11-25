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


def test_near_edges():
    """Grid points very close to an edge (but not on it) give accurate values.

    There, the weights of the two triangles that share the edge are large and
    nearly cancel. The value one float64 step from the midpoint of each edge
    between data points must agree with the mean of the values 1e-6 to
    either side of the edge.
    """
    rng = np.random.default_rng(5)
    x, y = rng.uniform(0, 100, (2, 10))
    z = rng.normal(size=(10, 1))
    grid = np.linspace(0, 100, 32)
    pts, vals, _ = ng._add_enclosing_triangle(x, y, z, grid, grid)

    def value(q):
        return ng._interpolate(pts, vals, q[:1], q[1:])[0, 0, 0]

    err = 0.0
    finite = True
    edges = set()
    for s in ng.Delaunay(pts).simplices:
        edges.update(tuple(sorted((s[a], s[b]))) for a, b in ((0, 1), (1, 2), (0, 2)))
    with np.errstate(divide='ignore', invalid='ignore'):
        for i, j in edges:
            if max(i, j) >= len(x):
                continue
            mid = (pts[i] + pts[j]) / 2
            normal = np.array([pts[i, 1] - pts[j, 1], pts[j, 0] - pts[i, 0]])
            normal *= 1e-6 / np.hypot(*normal)
            exp = (value(mid - normal) + value(mid + normal)) / 2
            for q in (mid, np.nextafter(mid, mid + normal)):
                v = value(q)
                finite &= bool(np.isfinite(v) and np.isfinite(exp))
                err = max(err, abs(v - exp))
    check('grid points close to edges', finite and err < 1e-7, '(max diff %.1e)' % err)


def test_edge_parallel_to_shift():
    """Grid points on an edge parallel to the first retry direction.

    Ill-conditioned grid points are retried at two points shifted along a
    fixed direction. If an edge is parallel to it, both shifted points are on
    the edge as well, and another direction must be used.
    """
    a = np.array([40.0, 40.0])
    b = a + 20 * np.array(ng.DIRECTIONS[0])
    pts = np.array([a, b, [20, 80], [80, 20], [10, 10], [90, 90],
                    [-1000, -1000], [-1000, 5000], [5000, -1000]])
    vals = np.array([[1.0], [3.0], [-2.0], [5.0], [0.5], [1.5], [0.0], [0.0], [0.0]])
    normal = np.array([a[1] - b[1], b[0] - a[0]])
    normal *= 1e-6 / np.hypot(*normal)

    def value(q):
        return ng._interpolate(pts, vals, q[:1], q[1:])[0, 0, 0]

    err = 0.0
    finite = True
    with np.errstate(divide='ignore', invalid='ignore'):
        for t in np.linspace(0.05, 0.95, 91):
            q = a + t * (b - a)
            v = value(q)
            exp = (value(q - normal) + value(q + normal)) / 2
            finite &= bool(np.isfinite(v) and np.isfinite(exp))
            err = max(err, abs(v - exp))
    check('edge parallel to the retry direction', finite and err < 1e-7, '(max diff %.1e)' % err)


def test_nonfinite_values():
    """Nonfinite data values are passed on, without retrying grid points."""
    rng = np.random.default_rng(6)
    x, y = rng.uniform(0, 100, (2, 10))
    z = rng.normal(size=(10, 2))
    z[3, 1] = np.nan
    grid = np.linspace(0, 100, 32)
    calls = []
    delaunay = ng.Delaunay
    ng.Delaunay = lambda p: calls.append(1) or delaunay(p)
    try:
        res = ng.natgrid(x, y, z, grid, grid)
    finally:
        ng.Delaunay = delaunay
    ok = (len(calls) == 1 and np.isfinite(res[:, :, 0]).all()
          and np.array_equal(res[:, :, 0], ng.natgrid(x, y, z[:, 0], grid, grid)))
    check('nonfinite values', ok, '(%d triangulations)' % len(calls))


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
    test_near_edges()
    test_edge_parallel_to_shift()
    test_nonfinite_values()
    test_noncontiguous()
    test_errors()
    if failures:
        print('%d test(s) FAILED' % len(failures))
        sys.exit(1)
    print('All tests passed')

"""Natural neighbor interpolation of scattered data to a rectangular grid.

This module replaces the C extension ``natgridmodule`` (NCAR's ``natgrid``,
based on Dave Watson's ``nngridr``), as far as PyBDSF uses it: natural
neighbor linear (Sibson) interpolation to a rectangular grid, with natgrid's
default parameters.

As in natgrid:

- The data points are triangulated together with the three vertices of a
  large enclosing triangle, which get the values of the least-squares plane
  through the data. They determine the interpolant outside (and close to) the
  convex hull of the data.
- An axis is rescaled if the extents of the region differ by more than a
  factor of two.
- The data coordinates are perturbed by a tiny pseudo-random amount (with the
  same random numbers as natgrid on Linux), to avoid degenerate cases.

The interpolant
---------------
The Sibson weight of a data point for a grid point q is the area that the
Voronoi cell of q (after inserting q) takes from the cell of the data point.
Watson's method collects these areas from the Delaunay triangles whose
circumcircle contains q. For such a triangle, with area A, barycentric
coordinates l0, l1, l2 of q, and D = R**2 - |q - c|**2 (c and R are the center
and radius of the circumcircle), the area given to vertex k has the closed
form ``D**2 * l_k / (16 * A * l0 * l1 * l2)``. Summed over the vertices, the
triangle adds::

    W = D**2 / (16 * A * l0 * l1 * l2)     to the total weight, and
    W * f(q)                               to the weighted sum of the values,

where f is the linear interpolant of the values over the triangle. The result
is the weighted sum divided by the total weight.
"""
import numpy as N
from scipy.spatial import Delaunay

# Number of grid rows that are processed at a time
BLOCK_ROWS = 32


def natgrid(x, y, z, xo, yo):
    """Interpolate scattered data to a rectangular grid.

    Parameters
    ----------
    x, y : array_like, shape (n,)
        Coordinates of the data points (at least four, all different).
    z : array_like, shape (n,) or (n, k)
        Data values; with shape (n, k), k sets of values are interpolated.
    xo, yo : array_like
        Increasing coordinates of the output grid along its two axes.

    Returns
    -------
    numpy.ndarray, shape (len(xo), len(yo)) or (len(xo), len(yo), k)
        Element [i, j] is the value at (xo[i], yo[j]).
    """
    x, y, xo, yo = (N.asarray(a, dtype=float).ravel() for a in (x, y, xo, yo))
    z = N.asarray(z, dtype=float)
    if not (x.size == y.size == len(z) >= 4) or z.ndim > 2:
        raise ValueError('x, y and z must have the same length (at least 4)')
    if xo.size == 0 or yo.size == 0 or N.any(N.diff(xo) < 0) or N.any(N.diff(yo) < 0):
        raise ValueError('The output coordinates must be increasing')
    zz = z.reshape(x.size, -1)

    pts, vals, scale = _add_enclosing_triangle(x, y, zz, xo, yo)
    with N.errstate(divide='ignore', invalid='ignore'):
        out = _interpolate(pts, vals, xo * scale[0], yo * scale[1])
    return out.reshape(xo.size, yo.size, *z.shape[1:])


def _add_enclosing_triangle(x, y, z, xo, yo):
    """Prepare the data points as natgrid does (``ReadDatad``).

    Returns the scaled and perturbed data points with the three vertices of
    the enclosing triangle appended, their values, and the scale factors.
    """
    # The region: the output grid, extended to include all data points.
    # (natgrid computes the limits in single precision.)
    def limits(v, vo):
        lo, hi, olo, ohi = (float(N.float32(a)) for a in (v.min(), v.max(), vo.min(), vo.max()))
        overlap = 1.01 * max(olo - lo, hi - ohi, 0.0)
        if overlap <= 1e-5:
            overlap = 0.01 * (ohi - olo)
        return olo - overlap, ohi + overlap

    (xmin, xmax), (ymin, ymax) = limits(x, xo), limits(y, yo)
    inside = (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)
    p, z = N.column_stack((x, y))[inside], z[inside]
    if len(p) < 4:
        raise ValueError('At least four data points near the output grid are needed')
    if len(N.unique(p, axis=0)) < len(p):
        raise ValueError('Two data points have the same coordinates')

    # Rescale the shorter axis if the extents differ by more than a factor 2
    lo = N.minimum(p.min(axis=0), (xmin, ymin))
    ext = N.maximum(p.max(axis=0), (xmax, ymax)) - lo
    scale = N.ones(2)
    ratio = ext[0] / ext[1]
    if ratio > 2 + 1e-5:
        scale[1] = ratio
    elif 1 / ratio > 2 + 1e-5:
        scale[0] = 1 / ratio
    p, lo, ext = p * scale, lo * scale, ext * scale

    # Enclosing triangle (60 times the extent), on the least-squares plane
    plane = N.linalg.lstsq(N.column_stack((N.ones(len(p)), p)), z, rcond=None)[0]
    big = lo + N.array([[-1, -1], [-1, 5], [5, -1]]) * ext * 10
    bigz = N.column_stack((N.ones(3), big)) @ plane

    # Tiny perturbation of the data coordinates
    r = _glibc_rand(367, p.size).reshape(p.shape) / 2147483647
    p = p + ext.max() * 1e-5 * (0.5 - r)
    return N.vstack((p, big)), N.vstack((z, bigz)), scale


def _glibc_rand(seed, count):
    """The first `count` values of glibc's rand() after srand(seed)."""
    r = [seed]
    for i in range(1, 31):
        r.append(16807 * r[-1] % 2147483647)
    r += r[:3]
    for i in range(34, 344 + count):
        r.append((r[i - 31] + r[i - 3]) & 0xFFFFFFFF)
    return N.array(r[344:], dtype=float) // 2


def _interpolate(pts, vals, xo, yo, fix=True):
    """Sibson interpolation of `vals` (shape (n, k)) at the grid (xo, yo)."""
    # Sorting the vertices by index makes the first one a data point (close
    # to the grid; the coordinates below are relative to it), and gives each
    # edge the same orientation in both triangles that share it.
    tri = N.sort(Delaunay(pts).simplices, axis=1)
    p0, p1, p2 = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
    z0, z1, z2 = vals[tri[:, 0]], vals[tri[:, 1]], vals[tri[:, 2]]

    # Circumcircles (center p0 + u, radius**2 r2) and twice the signed areas
    e1, e2 = p1 - p0, p2 - p0
    area2 = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]
    n1, n2 = (e1 ** 2).sum(axis=1), (e2 ** 2).sum(axis=1)
    u = N.column_stack((e2[:, 1] * n1 - e1[:, 1] * n2,
                        e1[:, 0] * n2 - e2[:, 0] * n1)) / (2 * area2[:, None])
    center = p0 + u
    r2 = (u ** 2).sum(axis=1)
    r = N.sqrt(r2)

    # Gradient of the linear interpolant f = z0 + gx * (x - x0) + gy * (y - y0)
    gx = ((z1 - z0) * e2[:, 1:] - (z2 - z0) * e1[:, 1:]) / area2[:, None]
    gy = ((z2 - z0) * e1[:, :1] - (z1 - z0) * e2[:, :1]) / area2[:, None]

    # The barycentric coordinates are l0 = L(p1, p2) / area2,
    # l1 = -L(p0, p2) / area2 and l2 = L(p0, p1) / area2, with L(a, b) the
    # line function of edge (a, b) below, so
    #     W = D**2 / (|area2| * l0 * l1 * l2) * constant
    #       = kappa * D**2 / (L(p1, p2) * L(p0, p2) * L(p0, p1)).
    # Near an edge, the weights of the two triangles that share it are large
    # and nearly cancel. L is computed in exactly the same way for both, so
    # that its rounding errors cancel as well.
    kappa = -area2 * N.abs(area2)
    # L(a, b) = (ay - by) * (x - ax) + (bx - ax) * (y - ay), as
    # (ay - by, bx - ax, ax, ay) for each triangle
    edges = [(a[:, 1] - b[:, 1], b[:, 0] - a[:, 0], a[:, 0], a[:, 1])
             for a, b in ((p1, p2), (p0, p2), (p0, p1))]

    def line(t, edge, xs, ys, out):
        cx, cy, ax, ay = (c[t] for c in edge)
        return N.add.outer(cx * (xs - ax), cy * (ys - ay), out=out)

    # Scratch arrays for one block, reused for every triangle (allocating
    # them for each triangle makes the interpolation about 10% slower)
    scratch = [N.empty(BLOCK_ROWS * yo.size) for _ in range(3)]
    total_buf = N.empty((BLOCK_ROWS, yo.size))
    wsum_buf = N.empty((BLOCK_ROWS, yo.size, vals.shape[1]))

    out = N.empty((xo.size, yo.size, vals.shape[1]))
    for r0 in range(0, xo.size, BLOCK_ROWS):
        xb = xo[r0:r0 + BLOCK_ROWS]
        total, wsum = total_buf[:xb.size], wsum_buf[:xb.size]
        total[:] = 0
        wsum[:] = 0
        for t in N.flatnonzero((center[:, 0] + r > xb[0]) & (center[:, 0] - r < xb[-1])):
            # The part of the block that the circumcircle covers
            i0, i1 = N.searchsorted(xb, center[t, 0] + (-r[t], r[t]))
            dxmin = max(xb[i0] - center[t, 0], center[t, 0] - xb[i1 - 1], 0) if i1 > i0 else r[t]
            h = N.sqrt(max(r2[t] - dxmin ** 2, 0))
            j0, j1 = N.searchsorted(yo, center[t, 1] + (-h, h))
            if i0 >= i1 or j0 >= j1:
                continue
            xs, ys = xb[i0:i1], yo[j0:j1]
            dx, dy = xs - p0[t, 0], ys - p0[t, 1]
            d, lines, f = (a[:xs.size * ys.size].reshape(xs.size, ys.size) for a in scratch)

            # w = kappa * D**2 / (product of the line functions), computed
            # in place in d (zero outside the circumcircle, where D <= 0)
            N.add.outer(dx * (2 * u[t, 0] - dx), dy * (2 * u[t, 1] - dy), out=d)
            N.maximum(d, 0, out=d)
            d *= d
            line(t, edges[0], xs, ys, out=lines)
            lines *= line(t, edges[1], xs, ys, out=f)
            lines *= line(t, edges[2], xs, ys, out=f)
            w = N.divide(d, lines, out=d, where=d > 0)
            w *= kappa[t]
            total[i0:i1, j0:j1] += w
            for k in range(vals.shape[1]):
                N.add.outer(gx[t, k] * dx + z0[t, k], gy[t, k] * dy, out=f)
                f *= w
                wsum[i0:i1, j0:j1, k] += f
        block = out[r0:r0 + BLOCK_ROWS]
        N.divide(wsum, total[:, :, None], out=block)

        # Grid points exactly on an edge (in practice never, because of the
        # perturbation) get the value at a slightly shifted point
        if fix and not N.isfinite(block).all():
            for i, j in N.argwhere(~N.isfinite(block).all(axis=2)):
                shift = 1e-9 * N.abs(pts).max()
                block[i, j] = _interpolate(pts, vals, xb[i:i + 1] + shift,
                                           yo[j:j + 1] + 0.755 * shift, fix=False)[0, 0]
    return out

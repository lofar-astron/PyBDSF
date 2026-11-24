"""Natural neighbor interpolation of scattered data to a rectangular grid.

This module replaces the C extension ``natgridmodule`` (from NCAR's ngmath
``natgrid`` package, based on Dave Watson's ``nngridr``), as far as PyBDSF
uses it: natural neighbor *linear* interpolation (Sibson interpolation) to a
rectangular grid, with the default control parameters of ``natgrid``
(no gradients, extrapolation allowed, automatic scaling, no duplicate input
points allowed).

Like ``natgrid``, the input points are triangulated together with the three
vertices of a large enclosing triangle, which get the values of the
least-squares plane through the data. These pseudo-points define the
interpolant outside the convex hull of the data (extrapolation) and close to
it. Also as in ``natgrid``, the input coordinates are perturbed by a tiny,
pseudo-random amount (the same amount as on Linux), and the axes are rescaled
if their extents differ by more than a factor of two.

How the interpolant is computed
-------------------------------
The Sibson coordinate of a data point for a grid point q is the area that the
Voronoi cell of q, after inserting q, takes from the Voronoi cell of the data
point. ``natgrid`` (Watson's method) obtains these areas from the Delaunay
triangles whose circumcircle contains q: for such a triangle T, with
circumcenter C, and for each vertex v of T, it computes the circumcenters g1
and g2 of q and the two edges of T that meet at v, and adds the signed area of
the triangle (g1, g2, C) to the weight of v.

The points C, g1 and g2 all lie on perpendicular bisectors of edges of T,
which gives a closed form for these areas. Let A be the area of T,
l0, l1 and l2 the barycentric coordinates of q with respect to T, and
D = R**2 - |q - C|**2, with R the circumradius of T (D > 0 means that q is
inside the circumcircle). Then the area added to the weight of vertex k is::

    D**2 * l_k / (16 * A * l0 * l1 * l2)

So each triangle T whose circumcircle contains q contributes
W_T = D**2 / (16 * A * l0 * l1 * l2) to the total weight, and
W_T * f_T(q) to the weighted sum of the data values, where f_T is the linear
interpolant of the data values over T (which is sum(l_k * z_k)). The
interpolated value is the ratio of these two sums. It only takes a few
operations per grid point and triangle, and these are done with vectorized
NumPy operations on blocks of the output grid, for the part of each block that
the circumcircle of a triangle covers.

The result is the same as that of ``natgrid`` up to rounding errors.
"""
import numpy as N
from scipy.spatial import Delaunay

# Constants of natgrid
_EPSILON = 0.00001
_RANGE = 10
_BIGTRI = ((-1, -1), (-1, 5), (5, -1))
_SEED = 367

# Number of grid rows (along the first output axis) processed at a time
_BLOCK_ROWS = 32


def natgrid(x, y, z, xo, yo):
    """Interpolate scattered data to a rectangular grid.

    Natural neighbor linear interpolation (Sibson interpolation), as done by
    ``natgridd`` of ``natgrid`` with its default parameters; see the module
    docstring.

    Parameters
    ----------
    x, y : array_like, shape (n,)
        Coordinates of the data points (at least four).
    z : array_like, shape (n,) or (n, k)
        Data values; with shape (n, k), k sets of values are interpolated at
        once (which is cheaper than interpolating them one by one).
    xo, yo : array_like, shape (nx,) and (ny,)
        Coordinates of the output grid along the first and second axis.
        Each must be monotonic (increasing or decreasing).

    Returns
    -------
    numpy.ndarray, shape (nx, ny) or (nx, ny, k)
        The interpolated values: element [i, j] is the value at
        (xo[i], yo[j]).

    Raises
    ------
    ValueError
        If the shapes do not match, if there are fewer than four data points
        (near the output grid), if an output axis is not monotonic, or if two
        data points have the same coordinates.
    """
    x = N.asarray(x, dtype=N.float64).ravel()
    y = N.asarray(y, dtype=N.float64).ravel()
    z = N.asarray(z, dtype=N.float64)
    xo = N.asarray(xo, dtype=N.float64).ravel()
    yo = N.asarray(yo, dtype=N.float64).ravel()
    if x.size != y.size or z.shape[0] != x.size or z.ndim > 2:
        raise ValueError('x, y and z must have the same number of data points')
    if x.size < 4:
        raise ValueError('At least four data points are needed')
    if xo.size == 0 or yo.size == 0:
        raise ValueError('The output grid must not be empty')
    single = z.ndim == 1
    z = z.reshape(x.size, -1)

    # The computation is done for increasing output coordinates
    xrev = xo.size > 1 and xo[0] > xo[-1]
    yrev = yo.size > 1 and yo[0] > yo[-1]
    if xrev:
        xo = xo[::-1]
    if yrev:
        yo = yo[::-1]
    if N.any(N.diff(xo) < 0) or N.any(N.diff(yo) < 0):
        raise ValueError('The output grid coordinates must be monotonic')

    pts, vals, magx, magy = _prepare(x, y, z, xo, yo)
    with N.errstate(divide='ignore', invalid='ignore'):
        out = _interpolate(pts, vals, xo * magx, yo * magy)

    if xrev:
        out = out[::-1]
    if yrev:
        out = out[:, ::-1]
    if single:
        out = out[:, :, 0]
    return N.ascontiguousarray(out)


def _prepare(x, y, z, xo, yo):
    """Select, scale and perturb the data points, and add the pseudo-points.

    This follows ``Initialized`` and ``ReadDatad`` of natgrid.

    Returns the (scaled) coordinates of the data points and the three
    pseudo-points (last), their values, and the scale factors of x and y.
    """
    # Limits of the output grid and of the data. natgrid computes these in
    # single precision (armind and armaxd).
    f32 = N.float32
    xstart, xend = float(f32(xo.min())), float(f32(xo.max()))
    ystart, yend = float(f32(yo.min())), float(f32(yo.max()))
    xil, xir = float(f32(x.min())), float(f32(x.max()))
    yib, yit = float(f32(y.min())), float(f32(y.max()))

    # Overlap regions: the smallest that include all data points
    horilap = vertlap = 0.0
    if xstart >= xil and xend <= xir:
        horilap = 1.01 * max(xstart - xil, xir - xend)
    elif xstart >= xil and xend >= xir:
        horilap = 1.01 * (xstart - xil)
    elif xstart <= xil and xend <= xir:
        horilap = 1.01 * (xir - xend)
    if horilap <= _EPSILON:
        horilap = 0.01 * (xend - xstart)
    if yib <= ystart and yend <= yit:
        vertlap = 1.01 * max(ystart - yib, yit - yend)
    elif ystart <= yib and yend <= yit:
        vertlap = 1.01 * (yit - yend)
    elif yib <= ystart and yit <= yend:
        vertlap = 1.01 * (ystart - yib)
    if vertlap <= _EPSILON:
        vertlap = 0.01 * (yend - ystart)
    minx, maxx = xstart - horilap, xend + horilap
    miny, maxy = ystart - vertlap, yend + vertlap

    # Only the data points inside the overlap region are used
    sel = (x > minx) & (x < maxx) & (y > miny) & (y < maxy)
    x, y, z = x[sel], y[sel], z[sel]
    n = x.size
    if n < 4:
        raise ValueError('At least four data points inside the overlap region '
                         'of the output grid are needed')

    # Extents of the data and overlap region; rescale an axis if they differ
    # by more than a factor of two
    magx = magy = 1.0
    for _ in range(2):
        lo = N.array([min((x * magx).min(), minx * magx),
                      min((y * magy).min(), miny * magy)])
        ext = N.array([max((x * magx).max(), maxx * magx),
                       max((y * magy).max(), maxy * magy)]) - lo
        qtxy = ext[0] / ext[1]
        qtyx = 1.0 / qtxy
        if qtxy > 2 + _EPSILON:
            magy *= qtxy
        elif qtyx > 2 + _EPSILON:
            magx *= qtyx
        else:
            break
    wbit = max(ext[0], ext[1]) * _EPSILON

    px = x * magx
    py = y * magy

    # Least-squares plane through the data (Cramer's rule, as in natgrid)
    sx, sy, sx2, sy2, sxy = px.sum(), py.sum(), (px * px).sum(), (py * py).sum(), (px * py).sum()
    sz = z.sum(axis=0)
    sxz = (px[:, None] * z).sum(axis=0)
    syz = (py[:, None] * z).sum(axis=0)
    det = (n * (sx2 * sy2 - sxy * sxy) - sx * (sx * sy2 - sy * sxy)
           + sy * (sx * sxy - sy * sx2))
    aaa = (sz * (sx2 * sy2 - sxy * sxy) - sxz * (sx * sy2 - sy * sxy)
           + syz * (sx * sxy - sy * sx2)) / det
    bbb = (n * (sxz * sy2 - syz * sxy) - sz * (sx * sy2 - sy * sxy)
           + sy * (sx * syz - sy * sxz)) / det
    ccc = (n * (sx2 * syz - sxy * sxz) - sx * (sx * syz - sy * sxz)
           + sz * (sx * sxy - sy * sx2)) / det

    # Vertices of the enclosing triangle, with the values of the plane
    big = lo + N.array(_BIGTRI, dtype=N.float64) * ext * _RANGE
    bigz = aaa + bbb * big[:, 0:1] + ccc * big[:, 1:2]

    # Duplicate data points are not allowed
    order = N.lexsort((py, px))
    same = (N.diff(px[order]) == 0) & (N.diff(py[order]) == 0)
    if N.any(same):
        i = N.flatnonzero(same)[0]
        raise ValueError('Data points %d and %d have the same coordinates'
                         % tuple(sorted(N.flatnonzero(sel)[order[i:i + 2]])))

    # Tiny perturbation of the coordinates, which natgrid uses to avoid
    # degenerate configurations
    r = _glibc_rand(_SEED, 2 * n).reshape(n, 2)
    px = px + wbit * (0.5 - r[:, 0] / 2147483647.0)
    py = py + wbit * (0.5 - r[:, 1] / 2147483647.0)

    pts = N.vstack((N.column_stack((px, py)), big))
    vals = N.vstack((z, bigz))
    return pts, vals, magx, magy


def _glibc_rand(seed, count):
    """Return the first `count` numbers of rand() of glibc after srand(seed).

    natgrid perturbs the data coordinates with rand(); this reproduces the
    sequence of glibc (the TYPE_3 additive feedback generator).
    """
    r = [0] * (344 + count)
    r[0] = seed
    for i in range(1, 31):
        r[i] = (16807 * r[i - 1]) % 2147483647
    for i in range(31, 34):
        r[i] = r[i - 31]
    for i in range(34, 344 + count):
        r[i] = (r[i - 31] + r[i - 3]) & 0xFFFFFFFF
    return N.array(r[344:], dtype=N.float64) // 2


def _interpolate(pts, vals, xo, yo, fix=True):
    """Sibson interpolation of `vals` at the grid points (xo[i], yo[j]).

    `pts` contains the data points followed by the vertices of an enclosing
    triangle, and `xo` and `yo` must be increasing.
    """
    # Sort the vertices of each triangle by index, so that the first one is
    # a data point (not a vertex of the enclosing triangle), and each edge
    # has the same first vertex in both triangles that share it
    tri = N.sort(Delaunay(pts).simplices, axis=1)
    p0, p1, p2 = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]

    # Twice the signed area, and the circumcircles. The coordinates are
    # relative to p0, which is close to the output grid.
    e1 = p1 - p0
    e2 = p2 - p0
    area2 = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]
    n1 = (e1 * e1).sum(axis=1)
    n2 = (e2 * e2).sum(axis=1)
    ux = (e2[:, 1] * n1 - e1[:, 1] * n2) / (2 * area2)
    uy = (e1[:, 0] * n2 - e2[:, 0] * n1) / (2 * area2)
    cx = p0[:, 0] + ux
    cy = p0[:, 1] + uy
    rad2 = ux * ux + uy * uy
    rad = N.sqrt(rad2)

    # The barycentric coordinate of vertex k is proportional to the line
    # function of the opposite edge (a, b):
    #     L(x, y) = (ay - by) * (x - ax) + (bx - ax) * (y - ay)
    # Close to an edge, the contributions of the two triangles that share it
    # are large and cancel. L is computed in the same way for both (with the
    # vertices in order of their index), so that its rounding errors cancel
    # as well.
    edges = []
    for a, b in ((p1, p2), (p0, p2), (p0, p1)):
        edges.append((a[:, 1] - b[:, 1], b[:, 0] - a[:, 0], a[:, 0], a[:, 1]))
    # l0 = L0 / area2, l1 = -L1 / area2, l2 = L2 / area2, so
    # W = D**2 / (|area2| * l0 * l1 * l2) = kappa * D**2 / (L0 * L1 * L2).
    # (W is the weight defined in the module docstring, times 8.)
    kappa = -area2 * N.abs(area2)

    # Linear interpolant over each triangle: f = z0 + fx * (x - x0) + fy * (y - y0)
    z0, z1, z2 = vals[tri[:, 0]], vals[tri[:, 1]], vals[tri[:, 2]]
    fx = ((z1 - z0) * e2[:, 1:2] - (z2 - z0) * e1[:, 1:2]) / area2[:, None]
    fy = ((z2 - z0) * e1[:, 0:1] - (z1 - z0) * e2[:, 0:1]) / area2[:, None]

    nx, ny = xo.size, yo.size
    nv = vals.shape[1]
    out = N.empty((nx, ny, nv), dtype=N.float64)
    den = N.empty((_BLOCK_ROWS, ny))
    num = N.empty((_BLOCK_ROWS, ny, nv))
    buf = [N.empty(_BLOCK_ROWS * ny) for _ in range(4)]

    for r0 in range(0, nx, _BLOCK_ROWS):
        r1 = min(r0 + _BLOCK_ROWS, nx)
        xb = xo[r0:r1]
        den[:] = 0.0
        num[:] = 0.0
        for t in N.flatnonzero((cx + rad > xb[0]) & (cx - rad < xb[-1])):
            # Rows and columns that the circumcircle of triangle t covers
            i0, i1 = N.searchsorted(xb, (cx[t] - rad[t], cx[t] + rad[t]))
            if i0 >= i1:
                continue
            dxmin = max(xb[i0] - cx[t], cx[t] - xb[i1 - 1], 0.0)
            h = N.sqrt(max(rad2[t] - dxmin * dxmin, 0.0))
            j0, j1 = N.searchsorted(yo, (cy[t] - h, cy[t] + h))
            if j0 >= j1:
                continue
            xs = xb[i0:i1]
            ys = yo[j0:j1]
            shape = (i1 - i0, j1 - j0)
            size = shape[0] * shape[1]
            d, l0, l1, l2 = (b[:size].reshape(shape) for b in buf)

            # D = R**2 - |q - c|**2 (zero outside the circumcircle), squared
            dx = xs - p0[t, 0]
            dy = ys - p0[t, 1]
            N.add.outer(dx * (2 * ux[t] - dx), dy * (2 * uy[t] - dy), out=d)
            N.maximum(d, 0.0, out=d)
            N.multiply(d, d, out=d)

            # Product of the line functions of the edges
            for lk, (bx, by, ax, ay) in zip((l0, l1, l2), edges):
                N.add.outer(bx[t] * (xs - ax[t]), by[t] * (ys - ay[t]), out=lk)
            N.multiply(l0, l1, out=l0)
            N.multiply(l0, l2, out=l0)

            # The weight W (zero outside the circumcircle)
            N.divide(d, l0, out=d, where=d > 0.0)
            d *= kappa[t]
            den[i0:i1, j0:j1] += d
            for v in range(nv):
                N.add.outer(fx[t, v] * dx + z0[t, v], fy[t, v] * dy, out=l1)
                l1 *= d
                num[i0:i1, j0:j1, v] += l1
        block = out[r0:r1]
        N.divide(num[:r1 - r0], den[:r1 - r0, :, None], out=block)
        if fix and not N.isfinite(block).all():
            for i, j in N.argwhere(~N.isfinite(block).all(axis=2)):
                block[i, j] = _shifted(pts, vals, xb[i], yo[j])
    return out


def _shifted(pts, vals, x, y):
    """Return the interpolated values at a point close to (x, y).

    At a grid point that lies exactly on an edge of the triangulation (or on
    the circumcircle of a triangle and on a line through one of its edges),
    the contributions of single triangles are infinite or undefined, although
    their sum is not. Such (very rare) points get the interpolated value at a
    slightly shifted position.
    """
    scale = N.abs(pts).max()
    for shift in (1e-9, 1e-7, 1e-5):
        # (in a direction that is unlikely to be along an edge)
        value = _interpolate(pts, vals, N.array([x + shift * scale]),
                             N.array([y + 0.7548776662 * shift * scale]), fix=False)[0, 0]
        if N.isfinite(value).all():
            return value
    return value

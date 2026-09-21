"""NumPy tessellation with the coordinate and label conventions of F2PY.

Pixel coordinates and generator labels start at one. Results are Fortran-order
float64 arrays. Spatial bounds only prune candidates; distances and ties use the
original Fortran arithmetic. The legacy sources are retained as test oracles.
"""

import operator

import numpy as np


def _inputs(n, m, xgens, ygens, snrgens, eps, code, ngens):
    n, m = operator.index(n), operator.index(m)
    if n <= 0 or m <= 0:
        raise ValueError("Image dimensions must be positive")
    x, y, snr = (np.asarray(a, dtype=np.float64) for a in (xgens, ygens, snrgens))
    if x.ndim != 1 or not x.size or y.shape != x.shape or snr.shape != x.shape:
        raise ValueError("Generator arrays must be nonempty, one-dimensional and equal length")
    if ngens is not None and operator.index(ngens) != x.size:
        raise ValueError("ngens must match the generator arrays")
    if not np.isfinite(x).all() or not np.isfinite(y).all() or not np.isfinite(eps):
        raise ValueError("Coordinates and eps must be finite")
    if isinstance(code, bytes):
        code = code.decode("ascii")
    if code not in ("s", "c"):
        raise ValueError("code must be 's' or 'c'")
    return n, m, x, y, float(eps), code


def _distance(px, py, x, y):
    dx, dy = px - x, py - y
    return np.sqrt(dx * dx + dy * dy)


def _bounds(lo_x, hi_x, lo_y, hi_y, x, y):
    near_x = np.maximum(np.maximum(lo_x - x, x - hi_x), 0.0)
    near_y = np.maximum(np.maximum(lo_y - y, y - hi_y), 0.0)
    far_x = np.maximum(np.abs(lo_x - x), np.abs(hi_x - x))
    far_y = np.maximum(np.abs(lo_y - y), np.abs(hi_y - y))
    return (np.sqrt(near_x * near_x + near_y * near_y),
            np.sqrt(far_x * far_x + far_y * far_y))


def _assign(n, m, x, y, weights, eps=0.0, fuzzy=False, centroids=None):
    if n * m <= 1024:
        return _assign_small(n, m, x, y, weights, eps, fuzzy, centroids)
    block = 16
    ni, nj = (n + block - 1) // block, (m + block - 1) // block
    image = np.empty((ni * block, nj * block), dtype=np.float64, order="F")
    tiles = image.reshape(ni, block, nj, block).transpose(0, 2, 1, 3)
    offsets = np.arange(block)
    for start in range(0, ni * nj, 512):
        ids = np.arange(start, min(start + 512, ni * nj))
        lo_x, lo_y = ids // nj * block + 1, ids % nj * block + 1
        hi_x, hi_y = np.minimum(lo_x + block - 1, n), np.minimum(lo_y + block - 1, m)
        lower, upper = _bounds(lo_x[:, None], hi_x[:, None], lo_y[:, None], hi_y[:, None], x, y)
        if centroids is None:
            lower /= weights
            upper /= weights
        else:
            cx, cy, inverse_radius = centroids
            near, far = _bounds(lo_x[:, None], hi_x[:, None], lo_y[:, None], hi_y[:, None], cx, cy)
            lower *= near * inverse_radius
            upper *= far * inverse_radius
            # A zero generator distance times an infinite radius is NaN at
            # singleton edge blocks. Such generators never win in Fortran.
            lower[:, ~np.isfinite(inverse_radius)] = np.inf
            upper[:, ~np.isfinite(inverse_radius)] = np.inf
        # A margin makes pruning conservative even at floating-point ties.
        cutoff = np.minimum(np.min(upper, axis=1), 1.e90)
        cutoff *= max(1.0, 1.0 + eps if fuzzy else 1.0)
        keep = lower <= (cutoff + abs(cutoff) * 1.e-12 + 1.e-300)[:, None]
        counts = keep.sum(axis=1)
        if np.any(counts == 0):
            raise ValueError("No generator has a defined distance below 1e90")
        # Batch tiles with equal candidate counts, preserving generator order.
        for count in np.unique(counts):
            selected = np.flatnonzero(counts == count)
            candidates = np.nonzero(keep[selected])[1].reshape(-1, count)
            tile_x, tile_y = ids[selected] // nj, ids[selected] % nj
            if count == 1 and np.all(upper[selected, candidates[:, 0]] < 1.e90):
                tiles[tile_x, tile_y] = 0 if fuzzy else candidates[:, 0, None, None] + 1
                continue
            px = np.minimum(lo_x[selected, None] + offsets, n)[:, :, None]
            py = np.minimum(lo_y[selected, None] + offsets, m)[:, None, :]
            nearest = np.full((len(selected), block, block), 1.e90)
            labels = np.zeros(nearest.shape, dtype=np.intp)
            for column in range(count):
                k = candidates[:, column, None, None]
                distance = _distance(px, py, x[k], y[k])
                if centroids is None:
                    distance /= weights[k]
                else:
                    # Keep both divisions: multiplication changes rounding at ties.
                    pixel_weights = 1.0 / (_distance(px, py, cx[k], cy[k]) * inverse_radius[k])
                    distance /= pixel_weights
                better = distance < nearest
                np.copyto(nearest, distance, where=better)
                np.copyto(labels, k + 1, where=better)
            if np.any(nearest >= 1.e90):
                raise ValueError("No generator has a defined distance below 1e90")
            if fuzzy:
                values = np.zeros(nearest.shape)
                threshold = (1.0 + eps) * nearest
                for column in range(count):
                    k = candidates[:, column, None, None]
                    distance = _distance(px, py, x[k], y[k]) / weights[k]
                    included = (distance <= threshold) & (labels != k + 1)
                    np.add(values, labels + k + 1, out=values, where=included)
                tiles[tile_x, tile_y] = values
            else:
                tiles[tile_x, tile_y] = labels
    return np.asfortranarray(image[:n, :m])


def _assign_small(n, m, x, y, weights, eps, fuzzy, centroids):
    """Avoid spatial bookkeeping for small images, using bounded row batches."""
    result = np.empty((n, m), dtype=np.float64, order="F")
    py = np.arange(1, m + 1, dtype=np.float64)[None, None, :]
    rows = max(1, 262144 // (m * len(x)))
    for start in range(0, n, rows):
        px = np.arange(start + 1, min(start + rows, n) + 1, dtype=np.float64)[None, :, None]
        distance = _distance(px, py, x[:, None, None], y[:, None, None])
        if centroids is None:
            distance /= weights[:, None, None]
        else:
            cx, cy, inverse_radius = (v[:, None, None] for v in centroids)
            pixel_weights = 1.0 / (_distance(px, py, cx, cy) * inverse_radius)
            distance /= pixel_weights
            distance[np.isnan(distance)] = np.inf
        winners = np.argmin(distance, axis=0)
        nearest = np.take_along_axis(distance, winners[None], axis=0)[0]
        if np.any(nearest >= 1.e90):
            raise ValueError("No generator has a defined distance below 1e90")
        if fuzzy:
            included = distance <= (1.0 + eps) * nearest
            np.put_along_axis(included, winners[None], False, axis=0)
            labels = np.arange(1, len(x) + 1)[:, None, None]
            result[start:start + rows] = np.sum(np.where(included, labels + winners + 1, 0), axis=0)
        else:
            result[start:start + rows] = winners + 1
    return result


def pytess_simple(n, m, xgens, ygens, snrgens, wts, eps, code, ngens=None):
    """Weighted Voronoi labels, or the legacy fuzzy sum when code is 'c'."""
    n, m, x, y, eps, code = _inputs(n, m, xgens, ygens, snrgens, eps, code, ngens)
    weights = np.asarray(wts, dtype=np.float64)
    if weights.shape != x.shape or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("Weights must be finite, positive and match the generators")
    if x.size == 1 and _bounds(1, n, 1, m, x, y)[1][0] / weights[0] < 1.e90:
        return np.full((n, m), 1.0 if code == "s" else 0.0, order="F")
    return _assign(n, m, x, y, weights, eps, code == "c")


def pytess_roundness(n, m, xgens, ygens, snrgens, eps, code, ngens=None):
    """The legacy two-pass roundness tessellation (hard labels only).

    Fuzzy roundness is rejected because the Fortran routine indexes arrays with
    fuzzy sums instead of generator labels, causing out-of-bounds accesses.
    """
    n, m, x, y, eps, code = _inputs(n, m, xgens, ygens, snrgens, eps, code, ngens)
    if code != "s":
        raise ValueError("Fuzzy roundness has no defined legacy result; use code='s'")
    if x.size == 1 and n * m > 1 and _bounds(1, n, 1, m, x, y)[1][0] < 1.e90:
        return np.ones((n, m), order="F")
    first = _assign(n, m, x, y, np.ones(x.size))
    labels = first.ravel(order="F").astype(np.intp) - 1
    with np.errstate(divide="ignore", invalid="ignore"):
        counts, cx, cy = _centroids(labels, n, m, x.size)
        dx = cx[labels].reshape(m, n)
        dy = cy[labels].reshape(m, n)
        dx -= np.arange(1, n + 1, dtype=np.float64)[None, :]
        dy -= np.arange(1, m + 1, dtype=np.float64)[:, None]
        np.square(dx, out=dx)
        np.square(dy, out=dy)
        dx += dy
        np.sqrt(dx, out=dx)
        # bincount accumulates in the same column-first order as Fortran.
        sums = np.bincount(labels, weights=dx.ravel(), minlength=x.size)
        inverse_radius = 1.0 / (sums / counts)
    valid = np.isfinite(inverse_radius)
    if not valid.any():
        raise ValueError("Roundness requires at least one tile with nonzero radius")
    # Empty and zero-radius tiles cannot win the reference's second pass.
    # Infinite bounds exclude them without evaluating NaN comparisons.
    inverse_radius[~valid] = np.inf
    cx[~valid] = cy[~valid] = -np.inf
    del first, labels, dx, dy
    with np.errstate(divide="ignore", invalid="ignore"):
        return _assign(n, m, x, y, None, centroids=(cx, cy, inverse_radius))


def _centroids(labels, n, m, count):
    if n * m * max(n, m) >= 2 ** 53:
        # Beyond exact integer sums, preserve the reference's accumulation order.
        area = np.bincount(labels, minlength=count)
        sx = np.bincount(labels, weights=np.tile(np.arange(1., n + 1), m), minlength=count)
        sy = np.bincount(labels, weights=np.repeat(np.arange(1., m + 1), n), minlength=count)
    else:
        # Each column contains runs of labels. Their integer coordinate sums can
        # be aggregated exactly, avoiding three floating reductions per pixel.
        change = np.empty(labels.size, dtype=bool)
        change[1:] = labels[1:] != labels[:-1]
        change[::n] = True
        starts = np.flatnonzero(change)
        lengths = np.diff(np.append(starts, labels.size))
        bins = labels[starts]
        sum_x = lengths * (2 * (starts % n + 1) + lengths - 1) // 2
        sum_y = lengths * (starts // n + 1)
        area = np.bincount(bins, weights=lengths, minlength=count)
        sx = np.bincount(bins, weights=sum_x, minlength=count)
        sy = np.bincount(bins, weights=sum_y, minlength=count)
    return area, sx / area, sy / area

"""NumPy tessellation with the coordinate and label conventions of F2PY.

Pixel coordinates and generator labels start at one. Results are Fortran-order
float64 arrays. Spatial bounds only prune candidates; distances and ties use the
original Fortran arithmetic. The legacy sources are retained as test oracles.
"""

import operator

import numpy as np


def _inputs(n, m, xgens, ygens, snrgens, eps, code, ngens):
    """Validate the shared F2PY-style arguments and normalize their types.

    Return ``(n, m, x, y, eps, code)`` with float64 coordinate vectors. SNR is
    checked for shape compatibility but does not affect either algorithm.
    ``ngens`` is optional for compatibility with the old generated wrapper;
    when supplied it must equal the coordinate-vector length.

    Invalid inputs raise ValueError rather than reproducing undefined Fortran
    behavior. Non-integral dimensions are rejected by ``operator.index``.
    """
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
    """Return Euclidean distances using broadcast-compatible coordinates.

    Keep the explicit squares, addition and square root: alternative formulas
    such as hypot or squared-distance comparisons can round differently at ties.
    """
    dx, dy = px - x, py - y
    return np.sqrt(dx * dx + dy * dy)


def _bounds(lo_x, hi_x, lo_y, hi_y, x, y):
    """Bound generator distances over closed, axis-aligned pixel rectangles.

    The lower bound is the distance to the nearest point of the rectangle
    (zero for a generator inside it). The upper bound uses the farthest corner.
    Inputs broadcast; block bounds shaped ``(blocks, 1)`` and generator vectors
    shaped ``(generators,)`` produce two ``(blocks, generators)`` arrays.
    """
    near_x = np.maximum(np.maximum(lo_x - x, x - hi_x), 0.0)
    near_y = np.maximum(np.maximum(lo_y - y, y - hi_y), 0.0)
    far_x = np.maximum(np.abs(lo_x - x), np.abs(hi_x - x))
    far_y = np.maximum(np.abs(lo_y - y), np.abs(hi_y - y))
    return (np.sqrt(near_x * near_x + near_y * near_y),
            np.sqrt(far_x * far_x + far_y * far_y))


def _assign(n, m, x, y, weights, eps=0.0, fuzzy=False, centroids=None):
    """Assign pixels after conservatively pruning generators for each block.

    With ``centroids=None``, use Euclidean distance divided by each positive
    generator weight. Otherwise ``centroids`` contains vectors of centroid x,
    centroid y and inverse mean radius; these define the legacy second-pass
    roundness score. The ``weights`` argument is then unused.

    For any block, the smallest candidate upper bound is an upper bound on the
    winning score everywhere. A candidate whose lower bound exceeds it cannot
    win. Fuzzy simple assignment expands that cutoff by ``1 + eps`` to retain
    all possible overlaps as well as the winner. Bounds decide what to skip;
    retained candidates still use the original floating-point expressions.

    Return a Fortran-contiguous float64 label/overlap map. Raise ValueError if
    any pixel has no usable score below the reference's initial 1e90 minimum.
    """
    if n * m <= 1024:
        return _assign_small(n, m, x, y, weights, eps, fuzzy, centroids)
    block = 16
    ni, nj = (n + block - 1) // block, (m + block - 1) // block
    image = np.empty((ni * block, nj * block), dtype=np.float64, order="F")
    # This is a view of the output, indexed as (block_x, block_y, pixel_x,
    # pixel_y). Batched assignments write directly into the final image.
    tiles = image.reshape(ni, block, nj, block).transpose(0, 2, 1, 3)
    offsets = np.arange(block)
    for start in range(0, ni * nj, 512):
        # Bound temporary arrays by processing at most 512 blocks at a time.
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
                # No pixel-level distances are needed. With no competitor,
                # the fuzzy sum is zero even though the hard label is nonzero.
                tiles[tile_x, tile_y] = 0 if fuzzy else candidates[:, 0, None, None] + 1
                continue
            # Repeat the last real coordinate in padding cells at image edges.
            # Pixel arrays broadcast to (selected_blocks, block, block).
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
                # Strict comparison keeps the earlier generator on exact ties.
                np.copyto(nearest, distance, where=better)
                np.copyto(labels, k + 1, where=better)
            if np.any(nearest >= 1.e90):
                raise ValueError("No generator has a defined distance below 1e90")
            if fuzzy:
                # Preserve the legacy encoding: add (winner + competitor)
                # for each overlap, excluding the winner itself. Labels are
                # one-based; a pixel without overlaps has value zero.
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
    """Evaluate all candidates directly when block bookkeeping would dominate.

    Distances have axes ``(generator, row, column)``. Row batching targets at
    most 262144 distances per batch, subject to a minimum of one row. Argmin
    selects the first generator on ties, matching the strict Fortran loop.
    The score definitions and return format are the same as in ``_assign``.
    """
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
            # Argmin would otherwise select a NaN; the reference's strict
            # comparison never accepts it as an improvement.
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
    """Compute weighted Voronoi labels or the legacy fuzzy overlap encoding.

    Parameters
    ----------
    n, m : int
        Positive image dimensions. Pixel coordinates are 1..n and 1..m.
    xgens, ygens : array-like, shape (ngens,)
        Finite generator coordinates in the same one-based coordinate system.
    snrgens : array-like, shape (ngens,)
        Unused SNR values, retained for compatibility with existing callers.
    wts : array-like, shape (ngens,)
        Positive finite divisors for Euclidean distances. Larger weights make
        generators more competitive; they are not squared-distance weights.
    eps : float
        Finite overlap allowance used only with code 'c'.
    code : {'s', 'c'}
        's' returns the nearest generator's one-based label. 'c' sums the
        winner's label plus each other label whose score is at most
        ``(1 + eps) * winner_score``. Without overlaps that sum is zero.
    ngens : int, optional
        Explicit generator count, which must match the input lengths.

    Returns
    -------
    numpy.ndarray, shape (n, m)
        Float64, Fortran-contiguous output, with integer-valued labels or sums.
        Exact ties choose the first generator in input order.
    """
    n, m, x, y, eps, code = _inputs(n, m, xgens, ygens, snrgens, eps, code, ngens)
    weights = np.asarray(wts, dtype=np.float64)
    if weights.shape != x.shape or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("Weights must be finite, positive and match the generators")
    if x.size == 1 and _bounds(1, n, 1, m, x, y)[1][0] / weights[0] < 1.e90:
        return np.full((n, m), 1.0 if code == "s" else 0.0, order="F")
    return _assign(n, m, x, y, weights, eps, code == "c")


def pytess_roundness(n, m, xgens, ygens, snrgens, eps, code, ngens=None):
    """The legacy two-pass roundness tessellation (hard labels only).

    Dimensions, coordinates, SNR and optional generator count follow
    ``pytess_simple``. ``eps`` is retained for signature compatibility but does
    not affect hard assignment; ``code`` must be 's'.

    First assign pixels with unit weights. For each resulting tile compute its
    centroid and inverse mean pixel distance to that centroid. Assign again
    with score ``distance_to_generator / (1 / (distance_to_centroid *
    inverse_mean_radius))``. The nested divisions preserve Fortran rounding.
    There are exactly two passes, not an iteration until convergence.

    Return a Fortran-contiguous float64 array of one-based labels. Empty and
    zero-radius first-pass tiles cannot win the second pass. Raise ValueError
    when no usable tile remains or a pixel has no valid score below 1e90.
    Fuzzy roundness is rejected because the Fortran routine indexes arrays with
    fuzzy sums instead of generator labels, causing out-of-bounds accesses.
    """
    n, m, x, y, eps, code = _inputs(n, m, xgens, ygens, snrgens, eps, code, ngens)
    if code != "s":
        raise ValueError("Fuzzy roundness has no defined legacy result; use code='s'")
    if x.size == 1 and n * m > 1 and _bounds(1, n, 1, m, x, y)[1][0] < 1.e90:
        return np.ones((n, m), order="F")
    first = _assign(n, m, x, y, np.ones(x.size))
    # Zero-based bin indices traverse x first, then y, exactly as the original
    # nested loops. Reshaping to (m, n) below preserves this traversal order.
    labels = first.ravel(order="F").astype(np.intp) - 1
    with np.errstate(divide="ignore", invalid="ignore"):
        counts, cx, cy = _centroids(labels, n, m, x.size)
        dx = cx[labels].reshape(m, n)
        dy = cy[labels].reshape(m, n)
        # Reuse the gathered centroid arrays as distance buffers to avoid
        # allocating full-image coordinate grids and squared-distance copies.
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
    """Return per-generator pixel counts and centroid x/y coordinates.

    ``labels`` is a flat, zero-based label map in Fortran traversal order;
    ``count`` includes generators with empty tiles. Empty centroids are NaN,
    which the caller handles when excluding unusable second-pass candidates.

    Integer pixel-coordinate sums below 2**53 are exactly representable in
    float64, so summing runs gives the same result as adding pixels individually.
    Above that conservative bound, fall back to the original addition order.
    """
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
        # Split at every column boundary too, so each run has constant y and
        # consecutive x values. Its x sum is an arithmetic progression.
        starts = np.flatnonzero(change)
        lengths = np.diff(np.append(starts, labels.size))
        bins = labels[starts]
        sum_x = lengths * (2 * (starts % n + 1) + lengths - 1) // 2
        sum_y = lengths * (starts // n + 1)
        area = np.bincount(bins, weights=lengths, minlength=count)
        sx = np.bincount(bins, weights=sum_x, minlength=count)
        sy = np.bincount(bins, weights=sum_y, minlength=count)
    return area, sx / area, sy / area

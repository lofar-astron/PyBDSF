"""Tessellation of the image into tiles, for the PSF variation.

This module replaces the Fortran extension `_pytesselate`, with the same
interface and the same results: every pixel gets the same value as in the
Fortran code, bit for bit.

Each pixel (i, j) belongs to the generator k with the smallest weighted
distance (in the Fortran code, pixel coordinates are 1-based, and ties go to
the generator with the lowest index):

    pytess_simple:     sqrt((i - x_k)**2 + (j - y_k)**2) / w_k
    pytess_roundness:  sqrt((i - x_k)**2 + (j - y_k)**2) / (1 / (r_k * c_k))

where, for `pytess_roundness`, r_k is the distance from the pixel to the
centroid of the tile of generator k in an unweighted tessellation, and c_k is
the inverse of the mean of r_k over that tile. If `code` is 's', the value of
a pixel is the (1-based) index of its tile. Otherwise ("fuzzy" tessellation),
it is the sum of (l + k) over all generators k != l that are within a factor
1 + eps of the distance to the nearest generator l (and 0 if there are none).

The Fortran code computes the distance to every generator for every pixel.
Here, the image is divided into blocks, and for each block, lower and upper
bounds of the distances are used to select the generators that can be the
nearest one (or within the fuzzy margin) for any of its pixels. For most
blocks, this is a single generator, which means that all of its pixels belong
to that tile. The other blocks are divided into smaller blocks, for which this
is repeated. Only for the pixels of small blocks with more than one candidate
are the distances computed, with exactly the same floating point operations
as in the Fortran code, so that the results are identical.
"""
import numpy as N

# Size of the large and small blocks, in pixels along each axis. _COARSE must
# be a multiple of _FINE.
_COARSE = 128
_FINE = 8

# If the number of pixels times the number of generators is at most this
# value, the whole image is processed as a single large block (which is
# faster for small images)
_DIRECT = 1 << 20

# Small blocks with different numbers of candidates are evaluated together
# (with padding) if that adds at most this number of distances
_PADDING = 4096

# Maximum number of distances computed at a time (number of small blocks
# times their number of candidates times the number of pixels per block)
_BATCH = 1 << 18

# Relative margin of the distance bounds, to make sure that rounding errors
# cannot exclude a generator that is (one of) the nearest for some pixel
_MARGIN = 1e-9

# Initial value of the minimum distance in the Fortran code: generators with
# a larger distance (or NaN) are never the nearest
_HUGE = 1e90

# Number of pixels processed at a time in _tile_properties
_CHUNK = 1 << 14


def pytess_simple(n, m, xgens, ygens, snrgens, wts, eps, code):
    """Weighted Voronoi tessellation of an n x m image.

    Parameters
    ----------
    n, m : int
        Shape of the image.
    xgens, ygens : array_like
        Positions of the generators, in pixels.
    snrgens : array_like
        Not used (kept for compatibility with the Fortran interface).
    wts : array_like
        Weights of the generators.
    eps : float
        Fractional overlap of the tiles for fuzzy tessellation.
    code : str
        's' for normal tessellation, anything else for fuzzy tessellation.

    Returns
    -------
    volrank : ndarray of float64, shape (n, m)
        The tile number (1-based) of each pixel, or for fuzzy tessellation,
        the sum of the pairs of tile numbers that a pixel belongs to.
    """
    xg, yg, wts = _as_arrays(xgens, ygens, wts)
    model = _Simple(xg, yg, wts, eps)
    return _tessellate(n, m, model, _is_fuzzy(code), eps)


def pytess_roundness(n, m, xgens, ygens, snrgens, eps, code):
    """Voronoi tessellation of an n x m image, modified to make the tiles rounder.

    First the image is tessellated without weights. Then it is tessellated
    again, with the distances to each generator scaled by the distance to the
    centroid of its tile, relative to the mean of that distance over the tile.

    The parameters and return value are the same as for `pytess_simple`
    (without the weights).

    In the Fortran code, fuzzy tessellation (code != 's') gave undefined
    results (or crashed), because the tile properties were computed from the
    fuzzy tessellation, in which pixel values are not tile numbers (and were
    used as array indices). Here they are computed from the normal
    tessellation.
    """
    xg, yg = _as_arrays(xgens, ygens)
    ngens = xg.size
    stats = tuple(N.zeros(ngens + 1) for _ in range(3))
    volrank = _tessellate(n, m, _Simple(xg, yg, N.ones(ngens), eps), False, eps,
                          stats=stats)
    xc, yc, roundpix = _tile_properties(volrank, stats)
    model = _Roundness(xg, yg, xc, yc, roundpix, eps)
    # The result of the first tessellation is not needed anymore, so its
    # memory is reused
    return _tessellate(n, m, model, _is_fuzzy(code), eps, out=volrank)


def _as_arrays(*arrays):
    arrays = [N.array(a, dtype=N.float64, ndmin=1).ravel() for a in arrays]
    if any(a.size != arrays[0].size for a in arrays):
        raise ValueError('all generator arrays must have the same length')
    return arrays


def _is_fuzzy(code):
    # The Fortran code uses only the first character of `code`
    return str(code)[:1] != 's'


def _pixel_dist(ii, jj, xg, yg):
    """sqrt((i - x)**2 + (j - y)**2) for the pixels of a set of blocks.

    `ii` and `jj` (shape (nblocks, npix)) are the pixel coordinates of the
    blocks, `xg` and `yg` (shape (nblocks, ncand)) the positions of the
    generators. Returns an array of shape (nblocks, ncand, npix, npix).
    """
    dx = ii[:, None, :] - xg[:, :, None]
    dx *= dx
    dy = jj[:, None, :] - yg[:, :, None]
    dy *= dy
    dist = dx[:, :, :, None] + dy[:, :, None, :]
    return N.sqrt(dist, out=dist)


class _Simple:
    """Distances for pytess_simple: sqrt((i - x)**2 + (j - y)**2) / w."""

    def __init__(self, xg, yg, wts, eps):
        self.xg, self.yg, self.wts = xg, yg, wts
        # The bounds are valid only for finite generators with positive
        # weights. Other generators (invalid input) are candidates for every
        # pixel, and if eps is not finite, all generators are.
        self.bounded = bool(N.isfinite(eps))
        self.regular = N.isfinite(xg) & N.isfinite(yg) & N.isfinite(wts) & (wts > 0)

    def candidates(self, ilo, ihi, jlo, jhi, fuzzy, eps, gens):
        """Generators that can matter for each block of a grid of blocks.

        The blocks are those of the rows `ilo` to `ihi` and the columns `jlo`
        to `jhi` (arrays of 1-based pixel coordinates of their first and
        last pixels), in row-major order. `gens` are the indices of the
        generators to consider. Returns a boolean array of shape (number of
        blocks, len(gens)).
        """
        cand = N.ones((ilo.size * jlo.size, gens.size), dtype=bool)
        if not self.bounded:
            return cand
        # If the nearest generator l is not regular, d_l is not larger than
        # the distance to any regular generator either, so it suffices to
        # select the regular generators among themselves
        regular = self.regular[gens]
        if N.all(regular):
            rgens = gens
        elif N.any(regular):
            rgens = gens[regular]
        else:
            return cand
        lo, hi = _dist_bounds(self.xg[rgens], self.yg[rgens], ilo, ihi, jlo, jhi)
        wts = self.wts[rgens]
        cand[:, regular] = _select(lo / wts, hi / wts, fuzzy, eps)
        return cand

    def distances(self, ii, jj, cand):
        """Weighted distances of the pixels of a set of blocks to their candidates.

        `ii` and `jj` (shape (nblocks, npix)) are the pixel coordinates of
        the blocks, and `cand` (shape (nblocks, ncand)) the indices of their
        candidates. Returns an array of shape (nblocks, ncand, npix, npix),
        and a function that gives the distance to the nearest generator l
        for the fuzzy comparison with each generator k, given the index of l
        in `cand` and the (weighted) distance to l.
        """
        dist = _pixel_dist(ii, jj, self.xg[cand], self.yg[cand])
        dist /= self.wts[cand][:, :, None, None]

        def distmin(nearest, best):
            return best[:, None]
        return dist, distmin


class _Roundness:
    """Distances for the second pass of pytess_roundness.

    In the Fortran code: dist = sqrt((i - x)**2 + (j - y)**2) / wts, where
    wts = 1 / (sqrt((xc - i)**2 + (yc - j)**2) * roundpix).
    """

    def __init__(self, xg, yg, xc, yc, roundpix, eps):
        self.xg, self.yg = xg, yg
        self.xc, self.yc, self.roundpix = xc, yc, roundpix
        self.bounded = bool(N.isfinite(eps))
        # Generators with an empty tile in the first pass have NaN distances,
        # and those with roundpix = inf (a tile of one pixel at its centroid)
        # have infinite or NaN distances. Neither can be the nearest
        # generator. For fuzzy tessellation, the latter belong to (almost)
        # every pixel, since inf <= (1 + eps) * inf. Generators that are not
        # finite (invalid input) are candidates for every pixel.
        finite = N.isfinite(xg) & N.isfinite(yg)
        self.regular = finite & N.isfinite(xc) & N.isfinite(yc) & N.isfinite(roundpix)
        self.always = ~finite
        self.always_fuzzy = ~finite | (N.isfinite(xc) & N.isfinite(yc) & (roundpix == N.inf))

    def candidates(self, ilo, ihi, jlo, jhi, fuzzy, eps, gens):
        """See _Simple.candidates."""
        if not self.bounded:
            return N.ones((ilo.size * jlo.size, gens.size), dtype=bool)
        cand = N.zeros((ilo.size * jlo.size, gens.size), dtype=bool)
        cand[:, self.always_fuzzy[gens] if fuzzy else self.always[gens]] = True
        regular = self.regular[gens]
        if not N.any(regular):
            return cand
        rgens = gens[regular]
        dlo, dhi = _dist_bounds(self.xg[rgens], self.yg[rgens], ilo, ihi, jlo, jhi)
        rlo, rhi = _dist_bounds(self.xc[rgens], self.yc[rgens], ilo, ihi, jlo, jhi)
        lo = dlo * rlo * self.roundpix[rgens]
        hi = dhi * rhi * self.roundpix[rgens]
        rcand = _select(lo, hi, False, eps)
        if fuzzy:
            # A generator k is added if d_k / wts_k <= (1 + eps) * d_l / wts_k,
            # where l is the nearest generator, and d the unweighted
            # distance. For finite, positive wts_k that requires
            # d_k <= (1 + eps) * d_l, approximately.
            dhi_max = N.where(rcand, dhi, 0.0).max(axis=1)
            thresh = max(1.0, 1.0 + eps) * dhi_max * (1.0 + _MARGIN)
            rcand |= dlo <= thresh[:, None]
        cand[:, regular] = rcand
        return cand

    def distances(self, ii, jj, cand):
        """See _Simple.distances."""
        dist0 = _pixel_dist(ii, jj, self.xg[cand], self.yg[cand])
        wts = _pixel_dist(ii, jj, self.xc[cand], self.yc[cand])
        wts *= self.roundpix[cand][:, :, None, None]
        N.divide(1.0, wts, out=wts)
        dist = dist0 / wts

        def distmin(nearest, best):
            # Fortran: the distance to the nearest generator, divided by the
            # weight of generator k (not of the nearest one)
            return N.take_along_axis(dist0, nearest[:, None], 1) / wts
        return dist, distmin


def _dist_bounds(xg, yg, ilo, ihi, jlo, jhi):
    """Lower and upper bounds of the distances of generators to a grid of blocks.

    See _Simple.candidates for the arguments. Returns arrays of shape
    (number of blocks, number of generators).
    """
    def bounds(lo, hi, g):
        # Per row or column of blocks: squared distance bounds along one
        # axis. With a = lo - g and b = g - hi, the minimum distance is
        # max(a, b, 0), and the maximum distance is |min(a, b)| (as hi >= lo).
        a = lo[:, None] - g
        b = g - hi[:, None]
        dmin = N.maximum(a, b)
        N.maximum(dmin, 0.0, out=dmin)
        dmax = N.minimum(a, b)
        return dmin * dmin, dmax * dmax
    dxlo, dxhi = bounds(ilo, ihi, xg)
    dylo, dyhi = bounds(jlo, jhi, yg)
    lo = dxlo[:, None, :] + dylo[None, :, :]
    hi = dxhi[:, None, :] + dyhi[None, :, :]
    shape = (-1, xg.size)
    return N.sqrt(lo, out=lo).reshape(shape), N.sqrt(hi, out=hi).reshape(shape)


def _select(lo, hi, fuzzy, eps):
    """Generators whose lower bound is within the margin of the minimum upper bound."""
    factor = max(1.0, 1.0 + eps) if fuzzy else 1.0
    thresh = hi.min(axis=1) * (factor * (1.0 + _MARGIN))
    return lo <= thresh[:, None]


def _fill_values(model, cand, gens, fuzzy):
    """Blocks that need exact evaluation, and the value of the other blocks.

    Blocks with a single (regular) candidate belong to that tile, and no other
    tile is within the fuzzy margin. Blocks without candidates occur only for
    invalid input, for which the Fortran result is undefined; they get 0.
    Returns the number of candidates, whether the block must be evaluated,
    and the value otherwise.
    """
    ncand = N.count_nonzero(cand, axis=1)
    first = gens[N.argmax(cand, axis=1)]
    single = (ncand == 1) & model.regular[first]
    mixed = (ncand > 1) | ((ncand == 1) & ~single)
    if fuzzy:
        fill = N.zeros(ncand.size)
    else:
        fill = N.where(single, first + 1.0, 0.0)
    return ncand, mixed, fill


def _tessellate(n, m, model, fuzzy, eps, out=None, stats=None):
    """Compute the tessellation of an n x m image with the given distance model.

    If `stats` is given (only for normal tessellation), it must be a tuple of
    three arrays of length ngens + 1, to which the number of pixels, and the
    sums of their i and j coordinates, are added for each tile number.
    """
    volrank = N.empty((n, m), order='F') if out is None else out
    gens = N.arange(model.xg.size)
    if gens.size == 0:
        # The Fortran result is undefined
        volrank[...] = 0.0
        return volrank
    if n * m * gens.size <= _DIRECT:
        # A single large block: its candidates are selected per small block
        with N.errstate(all='ignore'):
            result = _evaluate_block(model, 0, 0, _round_up(n), _round_up(m), gens, fuzzy, eps)
        volrank[...] = result[:n, :m]
        if stats is not None:
            _add_stats(stats, volrank, 0, 0)
        return volrank
    size = _COARSE
    a0 = N.arange(0, n, size)
    a1 = N.minimum(a0 + size, n)
    nblocks = a0.size
    # Pixel coordinates are 1-based in the Fortran code
    ilo = a0 + 1.0
    ihi = a1.astype(N.float64)
    with N.errstate(all='ignore'):
        for b0 in range(0, m, size):
            b1 = min(b0 + size, m)
            cand = model.candidates(ilo, ihi, N.array([b0 + 1.0]), N.array([float(b1)]),
                                    fuzzy, eps, gens)
            ncand, mixed, fill = _fill_values(model, cand, gens, fuzzy)
            # Runs of adjacent blocks with the same value are filled at once
            key = N.where(mixed, -1.0 - N.arange(nblocks), fill)
            starts = N.flatnonzero(N.r_[True, key[1:] != key[:-1]])
            ends = N.r_[starts[1:], nblocks]
            for k, e in zip(starts.tolist(), ends.tolist()):
                block = volrank[a0[k]:a1[e - 1], b0:b1]
                if mixed[k]:
                    result = _evaluate_block(model, a0[k], b0, size, size,
                                             N.flatnonzero(cand[k]), fuzzy, eps)
                    block[...] = result[:block.shape[0], :block.shape[1]]
                    if stats is not None:
                        _add_stats(stats, block, a0[k], b0)
                else:
                    block[...] = fill[k]
                    if stats is not None:
                        # Sums of integers, so they are exact
                        h, w = block.shape
                        t = int(fill[k])
                        stats[0][t] += h * w
                        stats[1][t] += w * (h * (2 * int(a0[k]) + h + 1) // 2)
                        stats[2][t] += h * (w * (2 * b0 + w + 1) // 2)
    return volrank


def _round_up(size):
    """`size` rounded up to a multiple of _FINE."""
    return -(-size // _FINE) * _FINE


def _evaluate_block(model, i0, j0, height, width, gens, fuzzy, eps):
    """Result for a large block of height x width pixels, with the given candidates.

    `i0` and `j0` are the (0-based) indices of its first pixel, and `height`
    and `width` must be multiples of _FINE. The block may extend beyond the
    image; those pixels are computed as well (and must be ignored by the
    caller).
    """
    nfi, nfj = height // _FINE, width // _FINE
    ilo = i0 + 1.0 + N.arange(0, height, _FINE)
    jlo = j0 + 1.0 + N.arange(0, width, _FINE)
    cand = model.candidates(ilo, ilo + (_FINE - 1), jlo, jlo + (_FINE - 1), fuzzy, eps, gens)
    ncand, mixed, fill = _fill_values(model, cand, gens, fuzzy)
    out = N.empty((nfi, _FINE, nfj, _FINE))
    out[...] = fill.reshape(nfi, 1, nfj, 1)
    mixed = N.flatnonzero(mixed)
    if mixed.size == 0:
        return out.reshape(height, width)
    pix = N.arange(_FINE)
    for group, count in _groups(mixed, ncand[mixed]):
        # Indices of the candidates of each block, in increasing order. Blocks
        # with fewer candidates are padded with repeats of their last one.
        counts = ncand[group]
        if counts[0] == count and N.all(counts == count):
            order = N.nonzero(cand[group])[1].reshape(group.size, count)
        else:
            order = N.argsort(~cand[group], axis=1, kind='stable')[:, :count]
            last = order[N.arange(group.size), counts - 1]
            order = N.where(N.arange(count) < counts[:, None], order, last[:, None])
        cgroup = gens[order]
        row, col = N.divmod(group, nfj)
        ii = ilo[row][:, None] + pix
        jj = jlo[col][:, None] + pix
        step = max(1, _BATCH // (count * _FINE * _FINE))
        for s in range(0, group.size, step):
            b = slice(s, s + step)
            out[row[b], :, col[b], :] = _evaluate(model, ii[b], jj[b], cgroup[b], fuzzy, eps)
    return out.reshape(height, width)


def _groups(blocks, ncand):
    """Divide blocks into groups, each with the maximum number of candidates in it.

    Blocks with the same number of candidates are in the same group, and
    groups are merged if the padding adds at most _PADDING distances.
    """
    counts = N.flatnonzero(N.bincount(ncand))[::-1].tolist()
    groups = []
    for count in counts:
        members = ncand == count
        if groups and N.count_nonzero(members) * (groups[-1][1] - count) * _FINE**2 <= _PADDING:
            groups[-1][0] |= members
        else:
            groups.append([members, count])
    return [(blocks[members], count) for members, count in groups]


def _evaluate(model, ii, jj, cand, fuzzy, eps):
    """Exact result for the pixels of a set of small blocks.

    `cand` (shape (nblocks, ncand)) are the indices of the candidates of each
    block, in increasing order, possibly followed by repeats of the last one.
    Returns an array of shape (nblocks, npix, npix).
    """
    dist, distmin = model.distances(ii, jj, cand)
    # Nearest generator, as in the Fortran code: the first one with a
    # distance smaller than all previous ones, and than _HUGE (so repeated
    # candidates do not matter)
    shape = dist.shape[:1] + dist.shape[2:]
    best = N.full(shape, _HUGE)
    tile = N.zeros(shape)  # 1-based index of the nearest generator, 0 if none
    nearest = N.zeros(shape, dtype=N.intp)  # its index in `cand`
    tiles = (cand + 1.0)[:, :, None, None]
    closer = N.empty(shape, dtype=bool)
    for k in range(cand.shape[1]):
        N.less(dist[:, k], best, out=closer)
        N.copyto(best, dist[:, k], where=closer)
        N.copyto(tile, tiles[:, k], where=closer)
        if fuzzy:
            N.copyto(nearest, k, where=closer)
    if not fuzzy:
        return tile
    # Fuzzy: sum of (l + k) for all k != l with dist_k <= (1 + eps) * dist_l
    within = dist <= (1.0 + eps) * distmin(nearest, best)
    # Repeated candidates must not be counted again
    within[:, 1:] &= (cand[:, 1:] != cand[:, :-1])[:, :, None, None]
    within &= nearest[:, None] != N.arange(cand.shape[1])[:, None, None]
    count = N.count_nonzero(within, axis=1)
    total = (tiles * within).sum(axis=1)
    # Where no generator was found, the Fortran result is undefined
    return N.where(tile > 0, count * tile + total, 0.0)


def _add_stats(stats, block, i0, j0):
    """Add the number of pixels of a block and their coordinate sums to `stats`."""
    npix, xsum, ysum = stats
    size = npix.size
    h, w = block.shape
    idx = block.astype(N.intp)
    # Number of pixels of each tile per row and per column. All sums are of
    # integers, so they are exact.
    rows = N.bincount((idx + size * N.arange(h)[:, None]).ravel(),
                      minlength=size * h).reshape(h, size)
    cols = N.bincount((idx + size * N.arange(w)).ravel(order='F'),
                      minlength=size * w).reshape(w, size)
    npix += rows.sum(axis=0)
    xsum += N.arange(i0 + 1.0, i0 + h + 1.0) @ rows
    ysum += N.arange(j0 + 1.0, j0 + w + 1.0) @ cols


def _tile_properties(volrank, stats):
    """Centroids of the tiles and the inverse of their mean distance to the centroid.

    This is the Fortran subroutine `tile_roundness`. `stats` holds the number
    of pixels of each tile and the sums of their coordinates (see
    `_tessellate`). The distances are summed in the same order as in the
    Fortran code, so the results are identical. Pixels with value 0 (which do
    not belong to any tile) are ignored.
    """
    n, m = volrank.shape
    npix, xsum, ysum = stats
    size = npix.size  # bin 0 is for pixels without a tile
    ii = N.arange(1, n + 1, dtype=N.float64)
    ncols = max(1, _CHUNK // n)
    with N.errstate(all='ignore'):
        xc = xsum / npix
        yc = ysum / npix
        # The distances must be summed sequentially, in the order of the
        # pixels in the Fortran code (i varies fastest). N.bincount does
        # that, so the sums are carried over from chunk to chunk by passing
        # them as the first weights.
        sumrad = N.zeros(size)
        bins = N.arange(size)
        for b0 in range(0, m, ncols):
            b1 = min(b0 + ncols, m)
            # The transpose of a Fortran-ordered block of columns is
            # C-contiguous, with the pixels in Fortran order
            idx = volrank[:, b0:b1].T.astype(N.intp)
            dx = xc[idx] - ii
            dx *= dx
            dy = yc[idx] - N.arange(b0 + 1.0, b1 + 1.0)[:, None]
            dy *= dy
            dx += dy
            N.sqrt(dx, out=dx)
            sumrad = N.bincount(N.concatenate((bins, idx.ravel())),
                                weights=N.concatenate((sumrad, dx.ravel())),
                                minlength=size)
        roundpix = 1.0 / (sumrad / npix)
    return xc[1:], yc[1:], roundpix[1:]

"""Tests for bdsf.tessellate, the replacement of the Fortran extension _pytesselate.

The results are compared, exactly, with a straightforward transcription of the
Fortran code (`reference_*` below), which computes the distance to every
generator for every pixel.

Run with: python ttessellate.py (exit status 0 if all tests pass)
"""
import sys
import warnings

import numpy as np

from bdsf import tessellate


# --- Reference implementation: a direct transcription of the Fortran code ---

def coords(n, m):
    """1-based pixel coordinates, with an axis for the generators."""
    i = np.arange(1, n + 1, dtype=float)[:, None, None]
    j = np.arange(1, m + 1, dtype=float)[None, :, None]
    return i, j


def reference_nearest(dist):
    """The nearest generator (1-based; 0 if none) as in the Fortran loop."""
    dumr = np.full(dist.shape[:2], 1e90)
    minind = np.zeros(dist.shape[:2], dtype=int)
    for k in range(dist.shape[2]):
        closer = dist[:, :, k] < dumr
        dumr[closer] = dist[:, :, k][closer]
        minind[closer] = k + 1
    return minind


def reference_fuzzy(dist, dist1, minind, eps):
    """Sum of (l + k) over the generators k != l within the fuzzy margin.

    dist1[..., k] is the distance to the nearest generator l, as used in the
    comparison with generator k.
    """
    volrank = np.zeros(minind.shape)
    for k in range(dist.shape[2]):
        add = (minind != k + 1) & (dist[:, :, k] <= (1.0 + eps) * dist1[:, :, k])
        volrank[add] += minind[add] + k + 1
    return volrank


def reference_simple(n, m, xg, yg, wts, eps, code):
    xg, yg, wts = (np.asarray(a, dtype=float) for a in (xg, yg, wts))
    i, j = coords(n, m)
    with np.errstate(all='ignore'):
        dist = np.sqrt((i - xg) * (i - xg) + (j - yg) * (j - yg)) / wts
    minind = reference_nearest(dist)
    if code == 's':
        return minind.astype(float)
    distmin = np.take_along_axis(dist, np.maximum(minind - 1, 0)[:, :, None], 2)
    return reference_fuzzy(dist, np.broadcast_to(distmin, dist.shape), minind, eps)


def reference_tile_roundness(minind, ngens):
    """Centroids and roundpix of the tiles; sums in the Fortran (column-major) order."""
    n, m = minind.shape
    npix = [0] * ngens
    x = [0.0] * ngens
    y = [0.0] * ngens
    for j in range(1, m + 1):
        for i in range(1, n + 1):
            ind = minind[i - 1, j - 1] - 1
            if ind >= 0:
                x[ind] += i
                y[ind] += j
                npix[ind] += 1
    with np.errstate(all='ignore'):
        x = [np.float64(x[k]) / npix[k] for k in range(ngens)]
        y = [np.float64(y[k]) / npix[k] for k in range(ngens)]
        sumrad = [0.0] * ngens
        for j in range(1, m + 1):
            for i in range(1, n + 1):
                ind = minind[i - 1, j - 1] - 1
                if ind >= 0:
                    sumrad[ind] += np.sqrt((x[ind] - i) * (x[ind] - i) +
                                           (y[ind] - j) * (y[ind] - j))
        roundpix = [1.0 / (np.float64(sumrad[k]) / npix[k]) for k in range(ngens)]
    return np.array(x), np.array(y), np.array(roundpix)


def reference_roundness(n, m, xg, yg, eps, code):
    xg, yg = (np.asarray(a, dtype=float) for a in (xg, yg))
    ngens = xg.size
    minind = reference_simple(n, m, xg, yg, np.ones(ngens), eps, 's').astype(int)
    x, y, roundpix = reference_tile_roundness(minind, ngens)
    i, j = coords(n, m)
    with np.errstate(all='ignore'):
        dumr1 = np.sqrt((x - i) * (x - i) + (y - j) * (y - j)) * roundpix
        wts = 1.0 / dumr1
        dist0 = np.sqrt((i - xg) * (i - xg) + (j - yg) * (j - yg))
        dist = dist0 / wts
        minind = reference_nearest(dist)
        if code == 's':
            return minind.astype(float)
        d_l = np.take_along_axis(dist0, np.maximum(minind - 1, 0)[:, :, None], 2)
        return reference_fuzzy(dist, d_l / wts, minind, eps)


# --- Helpers ---

def random_generators(rng, n, m, ngens, margin=0.0):
    xg = rng.uniform(-margin, n + margin, ngens)
    yg = rng.uniform(-margin, m + margin, ngens)
    return xg, yg


def check(result, expected, label=''):
    assert result.dtype == np.float64, label
    assert result.shape == expected.shape, label
    ndiff = np.count_nonzero(result != expected)
    assert ndiff == 0, f'{label}: {ndiff} of {result.size} pixels differ'


# (large block size, small block size, maximum padding)
BLOCK_SIZES = ((4, 2, 0), (12, 3, 1 << 30), (64, 8, 4096), 'default')


def for_block_sizes(test, sizes=BLOCK_SIZES):
    """Run `test` with several block sizes, so that many blocks have several
    candidates. Except for 'default', each image is divided into blocks, also
    if it is small, and the distances are computed for a few small blocks at
    a time. The maximum padding 0 means that small blocks are only evaluated
    together if they have the same number of candidates, 1 << 30 that they
    are all evaluated together.
    """
    names = ('_COARSE', '_FINE', '_PADDING', '_DIRECT', '_BATCH')
    old = [getattr(tessellate, name) for name in names]
    try:
        for size in sizes:
            if size == 'default':
                values = old
            else:
                values = list(size) + [0, 64]
            for name, value in zip(names, values):
                setattr(tessellate, name, value)
            test(size)
    finally:
        for name, value in zip(names, old):
            setattr(tessellate, name, value)


# --- Tests ---

def test_simple_random():
    rng = np.random.default_rng(1)

    def run(block):
        for n, m, ngens in ((70, 45, 1), (70, 45, 2), (64, 64, 7), (101, 37, 25)):
            xg, yg = random_generators(rng, n, m, ngens, margin=5.0)
            snr = rng.uniform(5, 500, ngens)
            for wts in (np.ones(ngens), np.log10(snr), np.sqrt(np.log10(snr))):
                for eps, code in ((0.05, 's'), (0.05, 'c'), (0.3, 'c'), (0.0, 'c')):
                    expected = reference_simple(n, m, xg, yg, wts, eps, code)
                    result = tessellate.tessellate_simple(n, m, xg, yg, snr, wts, eps, code)
                    check(result, expected, f'block {block}, {n}x{m}, {ngens} gens, {code}')
    for_block_sizes(run)


def test_simple_ties():
    """Generators at integer positions, many equal distances: lowest index wins."""
    def run(block):
        n = m = 41
        xg = [11.0, 31.0, 11.0, 31.0, 21.0, 21.0]
        yg = [11.0, 11.0, 31.0, 31.0, 21.0, 21.0]  # the last two coincide
        for wts in (np.ones(6), np.array([1.0, 1.0, 1.0, 1.0, 2.0, 2.0])):
            for code in 'sc':
                expected = reference_simple(n, m, xg, yg, wts, 0.05, code)
                result = tessellate.tessellate_simple(n, m, xg, yg, wts, wts, 0.05, code)
                check(result, expected, f'block {block}, {code}')
    for_block_sizes(run)
    result = tessellate.tessellate_simple(41, 41, [11.0, 31.0], [21.0, 21.0], [1, 1], [1, 1], 0.05, 's')
    assert np.all(result[:21] == 1) and np.all(result[21:] == 2)


def test_invalid_generators():
    """Weights and positions for which the bounds are not valid.

    These do not occur in PyBDSF. Such generators are candidates for every
    pixel, while the others are still selected by their bounds.
    """
    rng = np.random.default_rng(2)
    n, m, ngens = 60, 50, 12
    xg, yg = random_generators(rng, n, m, ngens)
    wts = rng.uniform(0.5, 2, ngens)
    cases = []
    for k, bad in ((3, 0.0), (3, -1.0), (3, np.inf), (3, np.nan), (0, -0.5)):
        w = wts.copy()
        w[k] = bad
        cases.append((xg, yg, w, 0.05))
    w = wts.copy()
    w[[1, 5, 9]] = (0.0, -2.0, np.inf)
    cases.append((xg, yg, w, 0.05))
    cases.append((xg, yg, -wts, 0.05))  # no valid weights at all
    for bad in (np.nan, np.inf):
        x = xg.copy()
        x[4] = bad
        cases.append((x, yg, wts, 0.05))
    cases.append((xg, yg, wts, np.nan))  # eps not finite

    def run(block):
        for x, y, w, eps in cases:
            i, j = coords(n, m)
            with np.errstate(all='ignore'):
                dist = np.sqrt((i - x) * (i - x) + (j - y) * (j - y)) / w
            # Where the Fortran code finds no generator, its result is undefined
            defined = reference_nearest(dist) > 0
            for code in 'sc':
                expected = reference_simple(n, m, x, y, w, eps, code)
                result = tessellate.tessellate_simple(n, m, x, y, w, w, eps, code)
                check(np.where(defined, result, 0.0), np.where(defined, expected, 0.0),
                      f'block {block}, w {w}, x {x}, eps {eps}, {code}')
        for bad in (np.nan, np.inf):
            x = xg.copy()
            x[4] = bad
            check(tessellate.tessellate_roundness(n, m, x, yg, x, 0.05, 's'),
                  reference_roundness(n, m, x, yg, 0.05, 's'), f'block {block}, x {x}')
    for_block_sizes(run)


def test_roundness_random():
    rng = np.random.default_rng(3)

    def run(block):
        for n, m, ngens in ((50, 35, 1), (50, 35, 3), (64, 48, 9), (37, 61, 20)):
            xg, yg = random_generators(rng, n, m, ngens)
            for eps, code in ((0.05, 's'), (0.05, 'c'), (0.2, 'c')):
                expected = reference_roundness(n, m, xg, yg, eps, code)
                result = tessellate.tessellate_roundness(n, m, xg, yg, xg, eps, code)
                check(result, expected, f'block {block}, {n}x{m}, {ngens} gens, {code}')
    for_block_sizes(run, sizes=((4, 2, 0), (64, 8, 1 << 30), 'default'))


def test_roundness_degenerate_tiles():
    """Empty tiles (NaN centroid) and tiles of one pixel (roundpix = inf)."""
    n, m = 30, 30
    cases = [
        # Duplicate generator: the second one has an empty tile
        ([5.0, 20.0, 20.0, 12.0], [5.0, 20.0, 20.0, 25.0]),
        # Generator outside the image, whose tile is empty
        ([5.0, 20.0, -100.0], [5.0, 20.0, 15.0]),
        # Generator in a corner, close to another one: a tile of one pixel
        ([1.0, 1.6, 20.0], [1.0, 1.6, 20.0]),
    ]

    def run(block):
        for xg, yg in cases:
            for code in 'sc':
                expected = reference_roundness(n, m, xg, yg, 0.05, code)
                result = tessellate.tessellate_roundness(n, m, xg, yg, xg, 0.05, code)
                check(result, expected, f'block {block}, {xg}, {code}')
    for_block_sizes(run, sizes=((4, 2, 1 << 30), 'default'))
    x, y, roundpix = tessellate._tile_properties(
        *_first_pass(n, m, *cases[0]))
    assert np.isnan(x[2]) and np.isnan(roundpix[2])
    x, y, roundpix = tessellate._tile_properties(
        *_first_pass(n, m, *cases[2]))
    assert roundpix[0] == np.inf


def _first_pass(n, m, xg, yg):
    xg, yg = np.asarray(xg), np.asarray(yg)
    stats = tuple(np.zeros(xg.size + 1) for _ in range(3))
    volrank = tessellate._tessellate(n, m, tessellate._Simple(xg, yg, np.ones(xg.size), 0.05),
                                     False, 0.05, stats=stats)
    return volrank, stats


def test_small_images():
    """Images smaller than a block, and of a single pixel or row."""
    rng = np.random.default_rng(4)

    def run(block):
        for n, m in ((1, 1), (1, 17), (17, 1), (5, 3), (33, 33)):
            xg, yg = random_generators(rng, n, m, 4)
            wts = rng.uniform(0.5, 2, 4)
            for code in 'sc':
                check(tessellate.tessellate_simple(n, m, xg, yg, wts, wts, 0.05, code),
                      reference_simple(n, m, xg, yg, wts, 0.05, code), f'{block}, {n}x{m}')
                check(tessellate.tessellate_roundness(n, m, xg, yg, wts, 0.05, code),
                      reference_roundness(n, m, xg, yg, 0.05, code), f'{block}, {n}x{m}')
    for_block_sizes(run)


def test_elongated_and_empty():
    """Very elongated images, and no generators (undefined in Fortran: 0 here)."""
    result = tessellate.tessellate_simple(1, 200000, [0.5], [7.0], [1.0], [1.0], 0.05, 's')
    assert result.shape == (1, 200000) and np.all(result == 1.0)
    xg, yg = [0.5, 1.5, 0.7], [10.0, 150000.0, 199990.0]
    check(tessellate.tessellate_simple(2, 200000, xg, yg, xg, [1.0, 2.0, 1.0], 0.05, 'c'),
          reference_simple(2, 200000, xg, yg, [1.0, 2.0, 1.0], 0.05, 'c'))
    check(tessellate.tessellate_roundness(2, 200000, xg, yg, xg, 0.05, 's'),
          reference_roundness(2, 200000, xg, yg, 0.05, 's'))
    for code in 'sc':
        assert np.all(tessellate.tessellate_simple(5, 4, [], [], [], [], 0.05, code) == 0)
        assert np.all(tessellate.tessellate_roundness(5, 4, [], [], [], 0.05, code) == 0)


def test_large_image():
    """A larger image, with the default block size."""
    rng = np.random.default_rng(5)
    n, m, ngens = 300, 260, 12
    xg, yg = random_generators(rng, n, m, ngens)
    wts = np.log10(rng.uniform(5, 500, ngens))
    for code in 'sc':
        check(tessellate.tessellate_simple(n, m, xg, yg, wts, wts, 0.05, code),
              reference_simple(n, m, xg, yg, wts, 0.05, code), code)
    check(tessellate.tessellate_roundness(n, m, xg, yg, wts, 0.05, 's'),
          reference_roundness(n, m, xg, yg, 0.05, 's'))


def test_interface():
    """Python lists and scalars are accepted, and only code[0] matters."""
    xg, yg, snr = [3.2, 10.7, 6.1], [4.4, 2.5, 11.9], [20.0, 30.0, 40.0]
    wts = [1.0, 1.5, 2.0]
    result = tessellate.tessellate_simple(13, 15, xg, yg, snr, wts, 0.05, 's')
    assert result.shape == (13, 15) and result.dtype == np.float64
    check(tessellate.tessellate_simple(13, 15, xg, yg, snr, wts, 0.05, 'simple'), result)
    fuzzy = tessellate.tessellate_simple(13, 15, xg, yg, snr, wts, 0.05, 'c')
    check(tessellate.tessellate_simple(13, 15, xg, yg, snr, wts, 0.05, 'x'), fuzzy)
    assert set(np.unique(result)) == {1.0, 2.0, 3.0}
    assert tessellate.tessellate_roundness(13, 15, xg, yg, snr, 0.05, 's').shape == (13, 15)
    try:
        tessellate.tessellate_simple(13, 15, xg, yg[:2], snr, wts, 0.05, 's')
    except ValueError:
        pass
    else:
        raise AssertionError('ValueError not raised for generators of different lengths')


def test_no_warnings():
    """No floating point warnings or errors, whatever the global NumPy settings."""
    old = np.seterr(all='raise')
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            tessellate.tessellate_simple(20, 20, [1.0, 1.0, 50.0], [1.0, 1.0, 50.0],
                                     [1, 1, 1], [1.0, 0.0, np.inf], 0.05, 'c')
            tessellate.tessellate_roundness(30, 30, [5.0, 20.0, 20.0, -100.0, 1.0, 1.6],
                                        [5.0, 20.0, 20.0, 15.0, 1.0, 1.6], [1] * 6, 0.05, 'c')
    finally:
        np.seterr(**old)


if __name__ == '__main__':
    failed = 0
    for name, test in list(globals().items()):
        if name.startswith('test_') and callable(test):
            try:
                test()
                print('PASSED', name)
            except Exception as e:
                failed += 1
                print('FAILED', name, repr(e))
    sys.exit(1 if failed else 0)

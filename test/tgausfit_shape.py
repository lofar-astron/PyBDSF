"""Unit checks for the rank test that decides whether a Gaussian shape can be fitted."""
import numpy as N

from bdsf.gausfit import _can_fit_gaussian_shape, _fit_beam


def mask_from_active(shape, pixels):
    """Return a mask (True = masked) with only `pixels` active."""
    mask = N.ones(shape, dtype=bool)
    for x, y in pixels:
        mask[x, y] = False
    return mask


def test_rank_check():
    shape = (12, 12)

    # A full block constrains every quadratic term.
    block = N.ones(shape, dtype=bool)
    block[2:6, 3:8] = False
    assert _can_fit_gaussian_shape(block)

    # Too few pixels.
    assert not _can_fit_gaussian_shape(mask_from_active(shape, [(1, 1), (2, 2), (3, 3)]))

    # A single row or a two-row strip cannot constrain the x curvature.
    assert not _can_fit_gaussian_shape(mask_from_active(shape, [(4, j) for j in range(10)]))
    assert not _can_fit_gaussian_shape(mask_from_active(shape, [(i, j) for i in (4, 5) for j in range(10)]))

    # A diagonal line has enough extent in x and y but is still degenerate.
    assert not _can_fit_gaussian_shape(mask_from_active(shape, [(i, i) for i in range(10)]))

    # A plus-shaped island: all points lie on the conic x*y = 0 (relative
    # to the centre), so the cross term is unconstrained.
    plus = [(5, j) for j in range(2, 9)] + [(i, 5) for i in range(2, 9)]
    assert not _can_fit_gaussian_shape(mask_from_active(shape, plus))

    # A 3x3 block is the smallest square that works.
    assert _can_fit_gaussian_shape(mask_from_active(shape, [(i, j) for i in range(3) for j in range(3)]))

    # Large islands are handled and well-conditioned.
    big = N.zeros((400, 400), dtype=bool)
    assert _can_fit_gaussian_shape(big)


def test_fit_beam():
    beam = (2.2, 2.0, 90.0)       # widened initial guess
    true_beam = (2.0, 2.0, 90.0)
    line = mask_from_active((10, 10), [(4, j) for j in range(10)])
    block = N.zeros((10, 10), dtype=bool)

    # Rank-deficient island: fit like fix_to_beam, with the true beam.
    assert _fit_beam(line, beam, true_beam, fix_to_beam=False) == (true_beam, True)
    # Shape can be fitted: free fit from the widened initial guess.
    assert _fit_beam(block, beam, true_beam, fix_to_beam=False) == (beam, False)
    # With fix_to_beam the existing behaviour is unchanged.
    assert _fit_beam(line, beam, true_beam, fix_to_beam=True) == (beam, True)


if __name__ == '__main__':
    test_rank_check()
    test_fit_beam()
    print('tgausfit_shape: all checks passed')

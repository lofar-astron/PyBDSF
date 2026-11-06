"""Tests for bdsf.mgfunction (MGFunction and lmder_fit)

Run with: python tmgfunction.py (exit status 0 if all tests pass)
"""
import sys
import warnings

import numpy as np

from bdsf.mgfunction import Gtype, MGFunction, lmder_fit


def gaussian(shape, a, x1, x2, s1, s2, pa):
    """Evaluate a 2D Gaussian (with the parametrisation of MGFunction)."""
    i, j = np.indices(shape, dtype=float)
    t = np.radians(pa)
    f1 = ((i - x1) * np.cos(t) + (j - x2) * np.sin(t)) / s1
    f2 = (-(i - x1) * np.sin(t) + (j - x2) * np.cos(t)) / s2
    return a * np.exp(-0.5 * (f1**2 + f2**2))


GAUSSIANS = [(5.0, 12.0, 10.0, 2.5, 1.7, 30.0), (3.0, 18.0, 15.0, 1.8, 1.6, 100.0)]
SHAPE = (30, 25)


def make_image(noise=0.0):
    image = sum(gaussian(SHAPE, *g) for g in GAUSSIANS)
    image += np.random.default_rng(1).normal(0, noise, SHAPE)
    return image.astype(np.float32)


def make_function(gtypes, params, image=None, mask=None):
    image = make_image(0.05) if image is None else image
    mask = (image < 0.3) if mask is None else mask
    fcn = MGFunction(image, mask, 1)
    for t, p in zip(gtypes, params):
        fcn.add_gaussian(t, p)
    return fcn, image, mask


def test_residual():
    params = [(4.0, 11.5, 10.3, 2.0, 2.2, 20.0), (2.5, 18.2, 14.8, 1.7, 1.6, 90.0)]
    fcn, image, mask = make_function([Gtype.g6, Gtype.g3], params)
    expected = (image - sum(gaussian(SHAPE, *p) for p in params))[~mask]
    residual = fcn.residual(fcn.get_fitted_parameters())
    assert np.allclose(residual, expected, rtol=0, atol=1e-12)
    assert np.isclose(fcn.chi2(), np.sum(expected**2), rtol=1e-12)


def test_jacobian():
    params = [(4.0, 11.5, 10.3, 2.0, 2.2, 20.0), (2.5, 18.2, 14.8, 1.7, 1.6, 90.0),
              (1.0, 5.0, 5.0, 1.5, 1.5, 0.0)]
    fcn, _, _ = make_function([Gtype.g6, Gtype.g3, Gtype.g1], params)
    x = fcn.get_fitted_parameters()
    assert len(x) == fcn.fitted_parameters() == 10
    jac = fcn.jacobian(x).copy()
    assert jac.shape == (10, fcn.data_size())
    for k in range(len(x)):
        h = 1e-6 * max(1.0, abs(x[k]))
        xp = x.copy()
        xp[k] += h
        xm = x.copy()
        xm[k] -= h
        numerical = (fcn.residual(xp).copy() - fcn.residual(xm)) / (2 * h)
        assert np.allclose(jac[k], numerical, rtol=0, atol=1e-8), k


def test_fit_recovers_gaussians():
    start = [(4.0, 11.5, 10.3, 2.0, 2.0, 20.0), (2.5, 18.2, 14.8, 1.7, 1.7, 90.0)]
    fcn, _, _ = make_function([Gtype.g6, Gtype.g6], start, image=make_image(),
                              mask=np.zeros(SHAPE, bool))
    assert lmder_fit(fcn, final=True, verbose=0)
    assert np.allclose(fcn.parameters, GAUSSIANS, rtol=1e-6)
    assert abs(fcn.find_peak()[0]) < 1e-6


def test_fixed_parameters():
    start = [(4.0, 11.5, 10.3, 2.0, 2.2, 20.0), (2.5, 18.2, 14.8, 1.7, 1.6, 90.0)]
    fcn, _, _ = make_function([Gtype.g3, Gtype.g1], start)
    assert lmder_fit(fcn, verbose=0)
    result = fcn.parameters
    assert result[0][3:] == start[0][3:]
    assert result[1][1:] == start[1][1:]
    assert result[0][:3] != start[0][:3]
    assert result[1][0] != start[1][0]


def test_repeated_fits():
    start = [(4.0, 11.5, 10.3, 2.0, 2.2, 20.0), (2.5, 18.2, 14.8, 1.7, 1.6, 90.0)]
    gtypes = [Gtype.g6, Gtype.g6]
    fcn, _, _ = make_function(gtypes, start)
    results = []
    for final in (False, True, False, True):
        fcn.reset()
        for t, p in zip(gtypes, start):
            fcn.add_gaussian(t, p)
        results.append((lmder_fit(fcn, final=final, verbose=0), fcn.parameters))
    # A new MGFunction does not have cached fit results
    for final in (False, True):
        fresh, _, _ = make_function(gtypes, start)
        results.append((lmder_fit(fresh, final=final, verbose=0), fresh.parameters))
    assert results[0] == results[2] == results[4]
    assert results[1] == results[3] == results[5]
    assert len(fcn._fit_results) == 2  # one for each value of `final`


def test_find_peak():
    image = make_image()
    mask = np.zeros(SHAPE, bool)
    mask[12, 10] = True  # mask the highest pixel
    fcn = MGFunction(image, mask)
    masked = np.where(mask, -np.inf, image)
    peak, (i, j) = fcn.find_peak()
    assert peak == masked.max()
    assert (i, j) == np.unravel_index(np.argmax(masked), SHAPE)
    fcn.add_gaussian(Gtype.g6, GAUSSIANS[0])
    peak, (i, j) = fcn.find_peak()
    assert (i, j) == (18, 15)


def test_python_interface():
    fcn = MGFunction(make_image(), np.zeros(SHAPE, bool))
    assert len(fcn) == 0 and fcn.fitted_parameters() == 0
    fcn.add_gaussian(Gtype.g6, GAUSSIANS[0])
    fcn.add_gaussian(Gtype.g3, [np.float32(1), 2, 3, 4, 5, 6])
    fcn.add_gaussian(Gtype.g1, GAUSSIANS[1])
    assert len(fcn) == 3 and fcn.fitted_parameters() == 10
    assert fcn[0] == GAUSSIANS[0] and fcn[-1] == GAUSSIANS[1]
    assert fcn[1] == (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    fcn[1] = (6, 5, 4, 3, 2, 1)
    assert fcn.parameters[1] == (6.0, 5.0, 4.0, 3.0, 2.0, 1.0)
    del fcn[1]
    assert fcn.parameters == GAUSSIANS and fcn.fitted_parameters() == 7
    fcn.remove_gaussian(-1)
    fcn.parameters = [GAUSSIANS[1]]
    assert fcn.parameters == [GAUSSIANS[1]] and fcn.errors == [(0.0,) * 6]
    fcn.reset()
    assert len(fcn) == 0 and fcn.parameters == []


def test_invalid_input():
    image = make_image()
    mask = np.zeros(SHAPE, bool)
    fcn = MGFunction(image, mask)
    for call, error in ((lambda: MGFunction(image[0], mask[0]), ValueError),
                        (lambda: MGFunction(image, mask[1:]), ValueError),
                        (lambda: MGFunction(image, mask.astype(int)), TypeError),
                        (lambda: MGFunction(image.astype(int), mask), TypeError),
                        (lambda: fcn.add_gaussian(Gtype.g6, (1, 2, 3)), ValueError),
                        (lambda: fcn.add_gaussian(2, GAUSSIANS[0]), ValueError),
                        (lambda: fcn[0], IndexError),
                        (lambda: fcn.remove_gaussian(0), IndexError)):
        try:
            call()
        except error:
            pass
        else:
            raise AssertionError('%s not raised' % error.__name__)


def test_degenerate_fits():
    # Fits that cannot work must fail without raising exceptions or issuing
    # warnings, even if numpy is set to raise on division by zero (as bdsf
    # does).
    old = np.seterr(divide='raise')
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            for params in ((1, 12, 10, 0.0, 1, 30), (1, 12, 10, 1, 1, np.inf),
                           (1, 12, 10, 1, 1, np.nan), (np.inf, 12, 10, 1, 1, 0)):
                fcn, _, _ = make_function([Gtype.g6], [params])
                assert not lmder_fit(fcn, verbose=0)
                fcn.find_peak()
                fcn.chi2()
            # Fewer pixels than fitted parameters
            mask = np.ones(SHAPE, bool)
            mask[12, 8:13] = False
            fcn, _, _ = make_function([Gtype.g6], [GAUSSIANS[0]], mask=mask)
            assert not lmder_fit(fcn, verbose=0)
            assert fcn.parameters == [GAUSSIANS[0]]
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

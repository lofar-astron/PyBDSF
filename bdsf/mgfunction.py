"""Multi-Gaussian function and Levenberg-Marquardt fitter.

Pure-Python (NumPy/SciPy) replacement of the former C++ extension module
``_cbdsm``.

An MGFunction with N 2D Gaussians is defined as::

    MGF(x1, x2) = sum_i A_i * exp(-0.5 * (f1_i**2 + f2_i**2))

    f1_i = ( (x1 - x1_i) * cos(pa_i) + (x2 - x2_i) * sin(pa_i)) / s1_i
    f2_i = (-(x1 - x1_i) * sin(pa_i) + (x2 - x2_i) * cos(pa_i)) / s2_i

Each Gaussian has the six parameters ``(A, x1, x2, s1, s2, pa)``, with the
position angle ``pa`` in degrees. Depending on its type (see `Gtype`), only
the first 1, 3 or all 6 of them are fitted; the others are kept fixed.

The function is only evaluated at the unmasked pixels of the data, with
vectorized NumPy operations. The results are cached, and only recomputed
when the parameters change.

The fits are ill-conditioned for some islands, so that tiny differences in
rounding can lead to different results. Therefore, the function and its
derivatives are evaluated with exactly the same floating point operations,
in the same order, as in the former C++ code (including the evaluation of the
exponential in long double precision), and the fits use MINPACK's LMDER, as
before (through SciPy). On x86-64, the results are therefore identical to
those of the C++ code, except for rare ill-conditioned fits: since version
1.15, SciPy uses a C translation of MINPACK, which in rare cases rounds
differently than the original Fortran code.
"""

from enum import IntEnum
from math import cos, pi, sin

import numpy as N
from scipy.optimize import leastsq

_DEG = pi / 180.0  # a Python float, like all per-Gaussian numbers


def _gaussian_scalars(s1, s2, pa):
    """Return the numbers, besides the parameters, that are needed to
    evaluate a Gaussian and its derivatives: (cos(pa), sin(pa),
    s1 / s2 - s2 / s1).

    For speed, these are computed with Python floats, which gives the same
    results as the former C++ code, except that Python raises an exception
    where C returns inf or nan: for the sine or cosine of infinity, and for a
    division by zero. These cases are handled separately.
    """
    pa *= _DEG
    try:
        cs = cos(pa)
        sn = sin(pa)
    except ValueError:  # pa is +/-inf
        cs = sn = N.nan
    try:
        ratio = s1 / s2 - s2 / s1
    except ZeroDivisionError:
        with N.errstate(all='ignore'):
            ratio = float(N.float64(s1) / s2 - s2 / N.float64(s1))
    return cs, sn, ratio


class Gtype(IntEnum):
    """Type of a Gaussian: the number of its fitted parameters."""
    g1 = 1  # amplitude only
    g3 = 3  # amplitude and position
    g6 = 6  # all parameters


class MGFunction:
    """Multi-Gaussian function, to be fitted with `lmder_fit`.

    Parameters
    ----------
    data : 2D ndarray of float32 or float64
        The image to fit.
    mask : 2D ndarray of bool
        Mask with the same shape as `data`; pixels where `mask` is True are
        ignored.
    weight : float
        Weight used in the chi^2 calculation.
    """

    def __init__(self, data, mask, weight=1.0):
        data = N.asarray(data)
        mask = N.asarray(mask)
        if data.ndim != 2 or mask.ndim != 2:
            raise ValueError("Data and mask should be rank-2 arrays")
        if data.shape != mask.shape:
            raise ValueError("Shape of mask doesn't match that of data")
        if mask.dtype != N.bool_:
            raise TypeError("Incorrect mask datatype")
        if data.dtype not in (N.float32, N.float64):
            raise TypeError("Incorrect data datatype")

        # Pixel coordinates and values of the unmasked pixels, in row-major
        # order.
        unmasked = ~mask
        self._x1, self._x2 = N.array(N.nonzero(unmasked), dtype=N.float64)
        self._data = data[unmasked].astype(N.float64)
        self._weight = float(weight)
        self._fit_results = {}  # see lmder_fit
        self.reset()

    def reset(self):
        """Forget all Gaussians."""
        self._gtypes = []
        self._params = N.empty((0, 6))
        self._layout_changed()

    def add_gaussian(self, gtype, parameters):
        """Add a Gaussian of type `gtype` (a `Gtype`) with the given six
        parameters."""
        gtype = Gtype(gtype)
        if len(parameters) != 6:
            raise ValueError("Wrong number of parameters for gaussian")
        self._gtypes.append(gtype)
        self._params = N.vstack((self._params, [float(p) for p in parameters]))
        self._layout_changed()

    @property
    def parameters(self):
        """Parameters of all Gaussians, as a list of 6-tuples."""
        return [tuple(p) for p in self._params.tolist()]

    def fitted_parameters(self):
        """Total number of fitted parameters."""
        return self._fit_idx.size

    def data_size(self):
        """Number of unmasked data points."""
        return self._data.size

    def find_peak(self):
        """Find the highest peak in the residual (data - function) image.

        Returns
        -------
        (value, (x1, x2))
        """
        with N.errstate(all='ignore'):
            resid = self._update()
        pidx = int(N.argmax(resid))
        if resid[pidx] != resid[pidx] and resid[0] == resid[0]:
            # Like the former C++ code, ignore NaNs unless the first value
            # is NaN
            pidx = int(N.nanargmax(resid))
        return float(resid[pidx]), (int(self._x1[pidx]), int(self._x2[pidx]))

    def chi2(self):
        """Chi^2 of the residual image (uniform weights)."""
        with N.errstate(all='ignore'):
            resid = self._update() / self._weight
            return float(N.dot(resid, resid))

    # ------------------------------------------------------------------
    # Interface for leastsq (see _lmder_fit). The fitted parameters are a
    # flat array, which excludes the fixed parameters of each Gaussian.
    # ------------------------------------------------------------------

    def get_fitted_parameters(self):
        """Return a copy of the fitted parameters."""
        return self._pflat[self._fit_idx]

    def set_fitted_parameters(self, x):
        """Set the fitted parameters."""
        self._pflat[self._fit_idx] = x

    def residual(self, x):
        """Return data - MGF for the fitted parameters `x`, for all unmasked
        pixels: an array of shape (data_size,), which must not be modified.
        """
        self.set_fitted_parameters(x)
        return self._update()

    def jacobian(self, x):
        """Return the derivatives of data - MGF with respect to the fitted
        parameters `x`: an array of shape (fitted_parameters, data_size),
        which must not be modified.
        """
        self.set_fitted_parameters(x)
        self._update()
        if self._jac is None:
            self._jac = self._compute_jacobian()
        return self._jac

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _layout_changed(self):
        """Update the bookkeeping after adding Gaussians or resetting.

        `_fit_idx` holds the indices of the fitted parameters (the first
        `gtype` parameters of each Gaussian) into `_pflat`, a flat view of
        the (ngaus, 6) parameter array.
        """
        self._fit_idx = N.array(
            [6 * i + j for i, t in enumerate(self._gtypes) for j in range(t)],
            dtype=N.intp)
        self._pflat = self._params.reshape(-1)
        self._cache_key = None

    def _update(self):
        """Evaluate the Gaussians for the current parameters, unless the
        cached results are for the same parameters already, and return
        data - MGF for all data points.

        This sets the following attributes:
        - _terms: for each Gaussian, a tuple (a, s1, s2, cs, sn, ratio, f1,
                  f2, val) with its parameters and the numbers returned by
                  _gaussian_scalars, and the values of f1, f2 and the
                  (unscaled) Gaussian for all data points
        - _resid: data - MGF for all data points
        - _jac:   None (the Jacobian is computed when needed)

        All values are computed exactly as in the former C++ code. The
        Gaussians are evaluated one by one, which is faster than evaluating
        them all at once for the small number of Gaussians in most islands.
        """
        key = self._params.tobytes()
        if key == self._cache_key:
            return self._resid
        resid = self._data.copy()
        terms = []
        for p in self._params.tolist():
            a, c1, c2, s1, s2, pa = p
            cs, sn, ratio = _gaussian_scalars(s1, s2, pa)
            dx1 = self._x1 - c1
            dx2 = self._x2 - c2
            f1 = (dx1 * cs + dx2 * sn) / s1
            f2 = (dx2 * cs - dx1 * sn) / s2
            # exp(-(f1**2 + f2**2) / 2), evaluated in long double precision
            val = (f1 * f1 + f2 * f2) * -0.5
            N.exp(val, out=val, dtype=N.longdouble)
            resid -= a * val
            terms.append((a, s1, s2, cs, sn, ratio, f1, f2, val))
        self._terms = terms
        self._resid = resid
        self._jac = None
        self._cache_key = key
        return resid

    def _compute_jacobian(self):
        """Return the derivatives of (data - MGF) with respect to the fitted
        parameters, as an array of shape (fitted_parameters, data_size).

        All values are computed exactly as in the former C++ code. The
        derivatives have a minus sign, because we differentiate data - MGF.
        """
        jac = N.empty((self._fit_idx.size, self._data.size))
        row = 0
        for gtype, (a, s1, s2, cs, sn, ratio, f1, f2, val) in zip(
                self._gtypes, self._terms):
            N.negative(val, out=jac[row])
            if gtype != Gtype.g1:
                v = val * -a
                # Position
                N.multiply(f1 * cs / s1 - f2 * sn / s2, v, out=jac[row + 1])
                N.multiply(f1 * sn / s1 + f2 * cs / s2, v, out=jac[row + 2])
                if gtype == Gtype.g6:
                    # Widths
                    N.divide(v * f1 * f1, s1, out=jac[row + 3])
                    N.divide(v * f2 * f2, s2, out=jac[row + 4])
                    # Position angle
                    N.multiply(v * _DEG * f1 * f2, ratio, out=jac[row + 5])
            row += gtype
        return jac


def lmder_fit(fcn, final=False, verbose=1):
    """Fit `fcn` to its data using the Levenberg-Marquardt algorithm.

    This uses LMDER from MINPACK (through scipy.optimize.leastsq), with an
    analytical Jacobian, like the former C++ code. On return, the parameters
    of `fcn` hold the result of the fit.

    Parameters
    ----------
    fcn : MGFunction
        The function to fit.
    final : bool
        If True, use a higher precision (for the final fit).
    verbose : bool
        If True, print a status line.

    Returns
    -------
    bool
        True if the fit converged.
    """
    final = bool(final)
    # The result of a fit only depends on the data (fixed for an MGFunction),
    # the types and initial parameters of the Gaussians, and `final`. When
    # gausfit fits an island, it usually repeats the same fit several times,
    # so the results are cached.
    key = (tuple(fcn._gtypes), fcn._params.tobytes(), final)
    if key not in fcn._fit_results:
        fcn._fit_results[key] = _lmder_fit(fcn, final)
    params, info, nfev, njev = fcn._fit_results[key]
    fcn._params[...] = params
    converged = 0 < info < 4

    if verbose:
        chi2 = fcn.chi2()
        print("status: %d  code: %d  Fev/Jev: %d/%d  chi2(/dp): %g(%g)  LMDER"
              % (converged, info, nfev, njev, chi2, chi2 / fcn.data_size()))

    return converged


def _lmder_fit(fcn, final):
    """Fit `fcn` with LMDER, and return (parameters, info, nfev, njev)."""
    npar = fcn.fitted_parameters()
    if npar == 0 or npar > fcn.data_size():
        # Improper input for LMDER
        return fcn._params.copy(), 0, 0, 0
    tol = 1e-6 if final else 1e-4
    x0 = fcn.get_fitted_parameters()
    # Like the former C++ code, ignore floating point errors (which may occur
    # when a fit diverges)
    with N.errstate(all='ignore'):
        x, _, infodict, _, info = leastsq(
            fcn.residual, x0, Dfun=fcn.jacobian, full_output=True,
            col_deriv=True, ftol=tol, xtol=tol, gtol=0.0, maxfev=200,
            factor=10.0)
    fcn.set_fitted_parameters(x)
    return fcn._params.copy(), info, infodict['nfev'], infodict['njev']

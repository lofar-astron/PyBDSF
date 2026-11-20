"""Multi-Gaussian function and Levenberg-Marquardt fitter.

Pure-Python (NumPy/SciPy) replacement of the former C++ extension module
``_cbdsm``.

An MGFunction with N 2D Gaussians is defined as::

    MGF(x1, x2) = sum_i amp_i * exp(-0.5 * (f1_i**2 + f2_i**2))

    f1_i = ( (x1 - x1_i) * cos_pa_i + (x2 - x2_i) * sin_pa_i) / sigma1_i
    f2_i = (-(x1 - x1_i) * sin_pa_i + (x2 - x2_i) * cos_pa_i) / sigma2_i

with, for Gaussian i:

- ``amp_i``: its peak value (amplitude)
- ``x1_i``, ``x2_i``: the pixel coordinates of its centre
- ``sigma1_i``, ``sigma2_i``: its widths (standard deviations, in pixels)
  along its two axes
- ``pa_i``: its position angle, in degrees, with ``cos_pa_i = cos(pa_i)``
  and ``sin_pa_i = sin(pa_i)``

so that each Gaussian has the six parameters
``(amp, x1, x2, sigma1, sigma2, pa)``. Depending on its type (see
`GaussianType`), only the first 1, 3 or all 6 of them are fitted; the others
are kept fixed.

The derivatives of a Gaussian ``g_i = amp_i * exp(-0.5 * (f1_i**2 +
f2_i**2))`` with respect to its parameters are::

    dg/d(amp)    = g / amp
    dg/d(x1)     = g * (f1 * cos_pa / sigma1 - f2 * sin_pa / sigma2)
    dg/d(x2)     = g * (f1 * sin_pa / sigma1 + f2 * cos_pa / sigma2)
    dg/d(sigma1) = g * f1**2 / sigma1
    dg/d(sigma2) = g * f2**2 / sigma2
    dg/d(pa)     = g * (pi / 180) * f1 * f2 * sigma12_term

with ``sigma12_term = sigma1 / sigma2 - sigma2 / sigma1``.

The function is only evaluated at the unmasked pixels of the data, with
vectorized NumPy operations. The results are cached, and only recomputed
when the parameters change.

The fits are ill-conditioned for some islands, so that tiny differences in
rounding can lead to different results. Therefore, the function and its
derivatives are evaluated with exactly the same floating point operations,
in the same order, as in the former C++ code (including the evaluation of the
exponential in long double precision), and the fits use the same algorithm
as before: scipy.optimize.leastsq, which uses MINPACK's LMDER. On x86-64, the
results are therefore identical to those of the C++ code, except for rare
ill-conditioned fits: since version 1.15, SciPy uses a C translation of
MINPACK, which in rare cases rounds differently than the original Fortran
code.
"""

from enum import IntEnum
from math import cos, pi, sin

import numpy as N
from scipy.optimize import leastsq

_DEG = pi / 180.0  # a Python float, like all per-Gaussian numbers


def _gaussian_scalars(sigma1, sigma2, pa):
    """Return the numbers, besides the parameters, that are needed to
    evaluate a Gaussian and its derivatives: (cos_pa, sin_pa, sigma12_term),
    see the module documentation.

    For speed, these are computed with Python floats, which gives the same
    results as the former C++ code, except that Python raises an exception
    where C returns inf or nan: for the sine or cosine of infinity, and for a
    division by zero. These cases are handled separately.
    """
    pa *= _DEG
    try:
        cos_pa = cos(pa)
        sin_pa = sin(pa)
    except ValueError:  # pa is +/-inf
        cos_pa = sin_pa = N.nan
    try:
        sigma12_term = sigma1 / sigma2 - sigma2 / sigma1
    except ZeroDivisionError:
        with N.errstate(all='ignore'):
            sigma12_term = float(N.float64(sigma1) / sigma2
                                 - sigma2 / N.float64(sigma1))
    return cos_pa, sin_pa, sigma12_term


class GaussianType(IntEnum):
    """Type of a Gaussian, which determines which of its parameters are
    fitted. Its value is the number of fitted parameters."""
    AMPLITUDE = 1  # amplitude only
    AMPLITUDE_AND_POSITION = 3  # amplitude and position (x1, x2)
    ALL_PARAMETERS = 6  # amplitude, position, widths and position angle


class MGFunction:
    """Multi-Gaussian function, to be fitted with `lm_fit_with_derivatives`.

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
        self._fit_results = {}  # see lm_fit_with_derivatives
        self.reset()

    def reset(self):
        """Forget all Gaussians."""
        self._gtypes = []
        self._params = N.empty((0, 6))
        self._layout_changed()

    def add_gaussian(self, gtype, parameters):
        """Add a Gaussian of type `gtype` (a `GaussianType`) with the given
        six parameters."""
        gtype = GaussianType(gtype)
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
            resid = self._evaluate_and_update_results()
        pidx = int(N.argmax(resid))
        if resid[pidx] != resid[pidx] and resid[0] == resid[0]:
            # Like the former C++ code, ignore NaNs unless the first value
            # is NaN
            pidx = int(N.nanargmax(resid))
        return float(resid[pidx]), (int(self._x1[pidx]), int(self._x2[pidx]))

    def residual_chi2(self):
        """Chi^2 of the residual (data - function) image, with uniform
        weights."""
        with N.errstate(all='ignore'):
            resid = self._evaluate_and_update_results() / self._weight
            return float(N.dot(resid, resid))

    # ------------------------------------------------------------------
    # Interface for leastsq (see _lm_fit_with_derivatives). The fitted
    # parameters are a flat array, which excludes the fixed parameters of
    # each Gaussian.
    # ------------------------------------------------------------------

    def get_fitted_parameters(self):
        """Return a copy of the fitted parameters."""
        return self._flat_params[self._fit_idx]

    def set_fitted_parameters(self, x):
        """Set the fitted parameters."""
        self._flat_params[self._fit_idx] = x

    def residual(self, x):
        """Return data - MGF for the fitted parameters `x`, for all unmasked
        pixels: an array of shape (data_size,), which must not be modified.
        """
        self.set_fitted_parameters(x)
        return self._evaluate_and_update_results()

    def jacobian(self, x):
        """Return the derivatives of data - MGF with respect to the fitted
        parameters `x`: an array of shape (fitted_parameters, data_size),
        which must not be modified.
        """
        self.set_fitted_parameters(x)
        self._evaluate_and_update_results()
        if self._jac is None:
            self._jac = self._compute_jacobian()
        return self._jac

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _layout_changed(self):
        """Update the bookkeeping after adding Gaussians or resetting.

        `_fit_idx` holds the indices of the fitted parameters (the first
        `gtype` parameters of each Gaussian) into `_flat_params`, a flat
        view of the (ngaus, 6) parameter array `_params`.
        """
        self._fit_idx = N.array(
            [6 * i + j for i, t in enumerate(self._gtypes) for j in range(t)],
            dtype=N.intp)
        self._flat_params = self._params.reshape(-1)
        self._cache_key = None

    def _evaluate_and_update_results(self):
        """Evaluate the Gaussians for the current parameters, unless the
        cached results are for the same parameters already, and return
        data - MGF for all data points.

        This sets the following attributes:
        - _terms: for each Gaussian, a tuple (amp, sigma1, sigma2, cos_pa,
                  sin_pa, sigma12_term, f1, f2, gaussian), with the
                  numbers defined in the module documentation, where f1, f2
                  and gaussian (the value of exp(-0.5 * (f1**2 + f2**2)),
                  without the amplitude) are arrays with their values for
                  all data points
        - _resid: data - MGF for all data points
        - _jac:   None (the Jacobian is computed when needed)

        All values are computed exactly as in the former C++ code. The
        Gaussians are evaluated one by one, which is faster than evaluating
        them all at once for the small number of Gaussians in most islands,
        and in place where possible, because temporary arrays slow down the
        fits of large islands.
        """
        key = self._params.tobytes()
        if key == self._cache_key:
            return self._resid
        resid = self._data.copy()
        terms = []
        for amp, x1_c, x2_c, sigma1, sigma2, pa in self._params.tolist():
            cos_pa, sin_pa, sigma12_term = _gaussian_scalars(sigma1, sigma2,
                                                             pa)
            dx1 = self._x1 - x1_c
            dx2 = self._x2 - x2_c
            # f1 = (dx1 * cos_pa + dx2 * sin_pa) / sigma1
            f1 = dx1 * cos_pa
            f1 += dx2 * sin_pa
            f1 /= sigma1
            # f2 = (-dx1 * sin_pa + dx2 * cos_pa) / sigma2
            f2 = dx2 * cos_pa
            f2 -= dx1 * sin_pa
            f2 /= sigma2
            # gaussian = exp((f1 * f1 + f2 * f2) / -2), with the exponential
            # evaluated in long double precision
            gaussian = f1 * f1
            gaussian += f2 * f2
            gaussian *= -0.5
            N.exp(gaussian, out=gaussian, dtype=N.longdouble)
            resid -= amp * gaussian
            terms.append((amp, sigma1, sigma2, cos_pa, sin_pa, sigma12_term,
                          f1, f2, gaussian))
        self._terms = terms
        self._resid = resid
        self._jac = None
        self._cache_key = key
        return resid

    def _compute_jacobian(self):
        """Return the derivatives of (data - MGF) with respect to the fitted
        parameters, as an array of shape (fitted_parameters, data_size).

        The derivatives are given in the module documentation. All values
        are computed exactly as in the former C++ code, in place (see
        _evaluate_and_update_results). The derivatives have a minus sign,
        because we differentiate data - MGF.
        """
        jac = N.empty((self._fit_idx.size, self._data.size))
        row = 0
        for gtype, terms in zip(self._gtypes, self._terms):
            (amp, sigma1, sigma2, cos_pa, sin_pa, sigma12_term,
             f1, f2, gaussian) = terms
            # Rows for this Gaussian (indexing jac creates a new view each time)
            d = list(jac[row:row + gtype])
            # Amplitude: -gaussian
            N.negative(gaussian, out=d[0])
            if gtype != GaussianType.AMPLITUDE:
                # -g, the value of the Gaussian with a minus sign
                neg_g = gaussian * -amp
                # Position: -g * (f1 * cos_pa / sigma1 - f2 * sin_pa / sigma2)
                #           -g * (f1 * sin_pa / sigma1 + f2 * cos_pa / sigma2)
                t = f1 * cos_pa
                t /= sigma1
                u = f2 * sin_pa
                u /= sigma2
                N.subtract(t, u, out=d[1])
                d[1] *= neg_g
                N.multiply(f1, sin_pa, out=t)
                t /= sigma1
                N.multiply(f2, cos_pa, out=u)
                u /= sigma2
                N.add(t, u, out=d[2])
                d[2] *= neg_g
                if gtype == GaussianType.ALL_PARAMETERS:
                    # Widths: -g * f1 * f1 / sigma1 and -g * f2 * f2 / sigma2
                    N.multiply(neg_g, f1, out=d[3])
                    d[3] *= f1
                    d[3] /= sigma1
                    N.multiply(neg_g, f2, out=d[4])
                    d[4] *= f2
                    d[4] /= sigma2
                    # Position angle: -g * deg * f1 * f2 * sigma12_term
                    N.multiply(neg_g, _DEG, out=d[5])
                    d[5] *= f1
                    d[5] *= f2
                    d[5] *= sigma12_term
            row += gtype
        return jac


def lm_fit_with_derivatives(function, final=False, verbose=1):
    """Fit `function` to its data using scipy.optimize.leastsq.

    This is the Levenberg-Marquardt algorithm, with an analytical Jacobian,
    like the former C++ code. On return, the parameters of `function` hold
    the result of the fit.

    Parameters
    ----------
    function : MGFunction
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
    key = (tuple(function._gtypes), function._params.tobytes(), final)
    if key not in function._fit_results:
        function._fit_results[key] = _lm_fit_with_derivatives(function, final)
    params, info, nfev, njev = function._fit_results[key]
    function._params[...] = params
    # Like the former C++ code, only treat info 1, 2 and 3 as convergence.
    # Because gtol is 0, info 4 (residual orthogonal to the Jacobian) only
    # occurs when the gradient is exactly 0, e.g. when all derivatives
    # underflow to 0 because the Gaussians moved away from the data.
    converged = 0 < info < 4

    if verbose:
        chi2 = function.residual_chi2()
        print("status: %d  code: %d  Fev/Jev: %d/%d  chi2(/dp): %g(%g)  LMDER"
              % (converged, info, nfev, njev, chi2,
                 chi2 / function.data_size()))

    return converged


def _lm_fit_with_derivatives(function, final):
    """Fit `function` with scipy.optimize.leastsq, and return (parameters,
    info, nfev, njev)."""
    npar = function.fitted_parameters()
    if npar == 0 or npar > function.data_size():
        # Improper input for leastsq
        return function._params.copy(), 0, 0, 0
    tol = 1e-6 if final else 1e-4
    x0 = function.get_fitted_parameters()
    # Like the former C++ code, ignore floating point errors (which may occur
    # when a fit diverges)
    with N.errstate(all='ignore'):
        x, _, infodict, _, info = leastsq(
            function.residual, x0, Dfun=function.jacobian, full_output=True,
            col_deriv=True, ftol=tol, xtol=tol, gtol=0.0, maxfev=200,
            factor=10.0)
    function.set_fitted_parameters(x)
    return function._params.copy(), info, infodict['nfev'], infodict['njev']

"""Multi-Gaussian function and Levenberg-Marquardt fitter.

Pure-Python (NumPy/SciPy) replacement of the former C++ extension module
``_cbdsm``.

An MGFunction with N 2D Gaussians is defined as::

    MGF(x1, x2) = sum_i A_i * exp(-0.5 * (f1_i**2 + f2_i**2))

    f1_i = ( (x1 - x1_i) * cos(pa_i) + (x2 - x2_i) * sin(pa_i)) / s1_i
    f2_i = (-(x1 - x1_i) * sin(pa_i) + (x2 - x2_i) * cos(pa_i)) / s2_i

Each Gaussian has the six parameters ``(A, x1, x2, s1, s2, pa)``: the
amplitude (``amp`` in the code), the centre (``center_x1``, ``center_x2``),
the widths (sigmas) along the major and minor axes (``sigma_major``,
``sigma_minor``), and the position angle ``pa`` of the major axis, in degrees.
``f1`` and ``f2`` are the offsets from the centre along the major and minor
axes, in units of sigma (``major_offset``, ``minor_offset``). The fit does
not enforce ``sigma_major >= sigma_minor``; gausfit takes care of that
afterwards. Depending on its type (see `GaussianType`), only the first 1, 3
or all 6 parameters are fitted; the others are kept fixed.

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


def _gaussian_scalars(sigma_major, sigma_minor, pa):
    """Return the numbers, besides the parameters, that are needed to
    evaluate a Gaussian and its derivatives: (cos(pa), sin(pa),
    sigma_major / sigma_minor - sigma_minor / sigma_major).

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
    # Only used in the derivative with respect to pa (see _compute_jacobian)
    try:
        sigma_ratio_diff = (sigma_major / sigma_minor
                            - sigma_minor / sigma_major)
    except ZeroDivisionError:
        with N.errstate(all='ignore'):
            sigma_ratio_diff = float(N.float64(sigma_major) / sigma_minor
                                     - sigma_minor / N.float64(sigma_major))
    return cos_pa, sin_pa, sigma_ratio_diff


class GaussianType(IntEnum):
    """Type of a Gaussian: the number of its fitted parameters."""
    amplitude_only = 1
    amplitude_and_position = 3
    all_parameters = 6


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
        """Add a Gaussian of type `gtype` (a `Gtype`) with the given six
        parameters."""
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
        """Chi^2 of the residual image (uniform weights)."""
        with N.errstate(all='ignore'):
            resid = self._evaluate_and_update_results() / self._weight
            return float(N.dot(resid, resid))

    # ------------------------------------------------------------------
    # Interface for leastsq (see _lm_fit_with_derivatives). The fitted
    # parameters are a flat array, which excludes the fixed parameters
    # of each Gaussian.
    # ------------------------------------------------------------------

    def get_fitted_parameters(self):
        """Return a copy of the fitted parameters."""
        return self._params_flat[self._fit_idx]

    def set_fitted_parameters(self, x):
        """Set the fitted parameters."""
        self._params_flat[self._fit_idx] = x

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
        `gtype` parameters of each Gaussian) into `_params_flat`, a flat
        view of the (ngaus, 6) parameter array.
        """
        self._fit_idx = N.array(
            [6 * i + j for i, t in enumerate(self._gtypes) for j in range(t)],
            dtype=N.intp)
        self._params_flat = self._params.reshape(-1)
        self._cache_key = None

    def _evaluate_and_update_results(self):
        """Evaluate the Gaussians for the current parameters, unless the
        cached results are for the same parameters already, and return
        data - MGF for all data points.

        This sets the following attributes:
        - _terms: for each Gaussian, a tuple (amp, sigma_major, sigma_minor,
                  cos_pa, sin_pa, sigma_ratio_diff, major_offset,
                  minor_offset, gaussian) with its parameters, the numbers
                  returned by _gaussian_scalars, and, for all data points,
                  the offsets from the centre along the major and minor axes
                  (in units of sigma) and the value of the (unscaled)
                  Gaussian
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
        for params in self._params.tolist():
            amp, center_x1, center_x2, sigma_major, sigma_minor, pa = params
            cos_pa, sin_pa, sigma_ratio_diff = _gaussian_scalars(
                sigma_major, sigma_minor, pa)
            dx1 = self._x1 - center_x1
            dx2 = self._x2 - center_x2
            # Offset along the major axis, in units of sigma:
            # (dx1 * cos_pa + dx2 * sin_pa) / sigma_major
            major_offset = dx1 * cos_pa
            major_offset += dx2 * sin_pa
            major_offset /= sigma_major
            # Offset along the minor axis, in units of sigma:
            # (-dx1 * sin_pa + dx2 * cos_pa) / sigma_minor
            minor_offset = dx2 * cos_pa
            minor_offset -= dx1 * sin_pa
            minor_offset /= sigma_minor
            # Unscaled Gaussian: exp((major_offset**2 + minor_offset**2) / -2),
            # with the exponential evaluated in long double precision
            gaussian = major_offset * major_offset
            gaussian += minor_offset * minor_offset
            gaussian *= -0.5
            N.exp(gaussian, out=gaussian, dtype=N.longdouble)
            resid -= amp * gaussian
            terms.append((amp, sigma_major, sigma_minor, cos_pa, sin_pa,
                          sigma_ratio_diff, major_offset, minor_offset,
                          gaussian))
        self._terms = terms
        self._resid = resid
        self._jac = None
        self._cache_key = key
        return resid

    def _compute_jacobian(self):
        """Return the derivatives of (data - MGF) with respect to the fitted
        parameters, as an array of shape (fitted_parameters, data_size).

        All values are computed exactly as in the former C++ code, in place
        (see _evaluate_and_update_results). The derivatives have a minus
        sign, because we differentiate data - MGF.
        """
        jac = N.empty((self._fit_idx.size, self._data.size))
        row = 0
        for gtype, (amp, sigma_major, sigma_minor, cos_pa, sin_pa,
                    sigma_ratio_diff, major_offset, minor_offset,
                    gaussian) in zip(self._gtypes, self._terms):
            # Rows of jac for the fitted parameters of this Gaussian, in the
            # order (amp, center_x1, center_x2, sigma_major, sigma_minor, pa)
            # (indexing jac creates a new view each time)
            rows = list(jac[row:row + gtype])
            N.negative(gaussian, out=rows[0])
            if gtype != GaussianType.amplitude_only:
                minus_amp_gaussian = gaussian * -amp
                # Centre, with g = minus_amp_gaussian:
                # x1: g * (major_offset * cos_pa / sigma_major
                #          - minor_offset * sin_pa / sigma_minor)
                # x2: g * (major_offset * sin_pa / sigma_major
                #          + minor_offset * cos_pa / sigma_minor)
                major_term = major_offset * cos_pa
                major_term /= sigma_major
                minor_term = minor_offset * sin_pa
                minor_term /= sigma_minor
                N.subtract(major_term, minor_term, out=rows[1])
                rows[1] *= minus_amp_gaussian
                N.multiply(major_offset, sin_pa, out=major_term)
                major_term /= sigma_major
                N.multiply(minor_offset, cos_pa, out=minor_term)
                minor_term /= sigma_minor
                N.add(major_term, minor_term, out=rows[2])
                rows[2] *= minus_amp_gaussian
                if gtype == GaussianType.all_parameters:
                    # Widths: g * major_offset**2 / sigma_major and
                    #         g * minor_offset**2 / sigma_minor
                    N.multiply(minus_amp_gaussian, major_offset, out=rows[3])
                    rows[3] *= major_offset
                    rows[3] /= sigma_major
                    N.multiply(minus_amp_gaussian, minor_offset, out=rows[4])
                    rows[4] *= minor_offset
                    rows[4] /= sigma_minor
                    # Position angle: g * deg * major_offset * minor_offset
                    #                 * sigma_ratio_diff
                    N.multiply(minus_amp_gaussian, _DEG, out=rows[5])
                    rows[5] *= major_offset
                    rows[5] *= minor_offset
                    rows[5] *= sigma_ratio_diff
            row += gtype
        return jac


def lm_fit_with_derivatives(func, final=False, verbose=1):
    """Fit `func` to its data using the Levenberg-Marquardt algorithm
    with derivatives (the Jacobian matrix).

    This uses scipy.optimize.leastsq, a wrapper around MINPACK’s lmdif and
    lmder algorithms. It uses an analytical Jacobian, like the former C++
    code. On return, the parameters of `func` hold the result of the fit.

    Parameters
    ----------
    func : MGFunction
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
    key = (tuple(func._gtypes), func._params.tobytes(), final)
    if key not in func._fit_results:
        func._fit_results[key] = _lm_fit_with_derivatives(func, final)
    params, info, nfev, njev = func._fit_results[key]
    func._params[...] = params
    converged = 0 < info < 4

    if verbose:
        chi2 = func.chi2()
        print("status: %d  code: %d  Fev/Jev: %d/%d  chi2(/dp): %g(%g)  LMDER"
              % (converged, info, nfev, njev, chi2, chi2 / func.data_size()))

    return converged


def _lm_fit_with_derivatives(func, final):
    """Fit `func` with LMDER, and return (parameters, info, nfev, njev)."""
    npar = func.fitted_parameters()
    if npar == 0 or npar > func.data_size():
        # Improper input for LMDER
        return func._params.copy(), 0, 0, 0
    tol = 1e-6 if final else 1e-4
    x0 = func.get_fitted_parameters()
    # Like the former C++ code, ignore floating point errors (which may occur
    # when a fit diverges)
    with N.errstate(all='ignore'):
        x, _, infodict, _, info = leastsq(
            func.residual, x0, Dfun=func.jacobian, full_output=True,
            col_deriv=True, ftol=tol, xtol=tol, gtol=0.0, maxfev=200,
            factor=10.0)
    func.set_fitted_parameters(x)
    return func._params.copy(), info, infodict['nfev'], infodict['njev']

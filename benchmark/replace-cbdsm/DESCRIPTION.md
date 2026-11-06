# Replace the C++ extension `_cbdsm` with Python

Branch: `replace_cbdsm`, on top of `master` (4fc79a2).

## What changed

- **New module `bdsf/mgfunction.py`** with `MGFunction`, `Gtype` and
  `lmder_fit`, which have the same interface as those of `_cbdsm` (the
  unused `dn2g_fit` and `dnsg_fit` were dropped). `gausfit.py` imports them
  from there; nothing else in `gausfit.py` changed.
- **Fitting:** `lmder_fit` uses MINPACK's LMDER through
  `scipy.optimize.leastsq` (the same algorithm, with the same settings, as
  the C++ code), with an analytical Jacobian.
- **Exact arithmetic:** the fits of some islands are ill-conditioned: a 1-ulp
  difference in the function or its Jacobian can lead to a different
  result. A first, faster version (one matrix multiplication for all
  Gaussians, float64 `exp`) changed 2 Gaussians of the test image by > 3%,
  and the test failed. The Gaussians and their derivatives are therefore
  computed with exactly the same floating point operations, in the same
  order, as in the C++ code: Python's `math.cos`/`math.sin` (same libm),
  true divisions, no fused multiply-adds, and the exponential in long double
  precision (`N.exp(x, dtype=N.longdouble, out=float64_array)`). Replaying
  recorded fits (see below) with the original Fortran LMDER (built with f2py)
  and this evaluation reproduces the C++ results bit for bit in all cases
  that were checked.
- **Speed:** the Gaussians are evaluated one by one with NumPy (over all
  pixels), which has the least overhead for the 1–3 Gaussians of most
  islands. Values are cached per parameter set, so the Jacobian reuses the
  Gaussians of the preceding function evaluation (and the extra evaluations
  that `leastsq` does before calling LMDER cost almost nothing; calling
  SciPy's private LMDER wrapper directly was only 2% faster). Still, a single fit of a
  small island is several times slower than in C++ (NumPy has ~0.5 µs
  overhead per operation, and an evaluation needs ~40 operations).
  This is more than made up for by caching the **results** of fits:
  `fit_island` repeats `fit_iter` up to 15 times, and most fits are exact
  repeats of an earlier one (75% of the fits for the test image, 84% for
  NGC6048), because their input did not change. `lmder_fit` caches the
  results per `MGFunction` (keyed on the Gaussian types, the initial
  parameters and `final`), which gives exactly the same results.
- **Floating point errors** in diverging fits (division by zero, overflow)
  give inf/nan as in C++, regardless of the global NumPy error settings
  (`gaul2srl.py` sets `divide='raise'`).
- **Build:** `src/c++`, `src/minpack` and `src/port3` were removed, as were
  Boost (CMake, `cibuildwheel/before_*.sh`, the cibuildwheel environment in
  `pyproject.toml`, the CI workflow, README) and C++ (the CMake project
  languages). Fortran (`_pytesselate`) and C (`natgrid`) remain.
- **New test** `test/tmgfunction.py` (added to CI and the cibuildwheel
  `test-command`): residuals and Jacobian (against finite differences) for
  mixed Gaussian types, fitting, fixed parameters, cached fits, `find_peak`,
  the Python interface, invalid input, and degenerate fits (no exceptions or
  warnings, also with `seterr(divide='raise')`).

## Results

Machine: 20-core i7-12700, NumPy 2.4.6, SciPy 1.18.1. Baseline: the
`_cbdsm` from a wheel build of `master`.

**Single fits** (`replay_fits.py`, best of 3, no result cache):

| Recorded fits | Bit-identical | C++ | Python | Ratio |
| --- | --- | --- | --- | --- |
| test image (831) | 831 | 0.48 s | 1.72 s | 3.6 |
| NGC6048 (32983) | 32979 | 39.4 s | 45.9 s | 1.16 |

For NGC6048, the ratio is 7.2 for islands with < 50 pixels, and 0.91 for
islands with ≥ 1000 pixels.

**End-to-end** (`run_comparison.sh`, 3 interleaved runs each):

| | C++ | Python |
| --- | --- | --- |
| `test/tbdsf_process_image.py` (2 cores) | 5.70–5.78 s | 5.27–5.36 s |
| `run_pybdsf.py` NGC6048 (20 cores), total | 64.0–64.1 s (one run: 74.1 s) | 54.9–55.3 s |
| … of which `Op_gausfit` | 12.8–12.9 s (one run: 22.8 s) | 7.7–8.0 s |
| … of which `Op_wavelet_atrous` | 40.2–40.4 s | 36.3–36.5 s |

- Test image: all 50 FITS outputs (images and catalogs) are identical to
  those of C++.
- NGC6048: all 1077 Gaussians of the image are identical, as are the number
  of islands and Gaussians on each wavelet scale. 2 of the 415 wavelet
  Gaussians differ by > 1e-6 (by 47% and 33%); the others differ by < 1e-13
  or not at all.

## Remaining difference

The 4 non-identical fits of NGC6048 (all ill-conditioned; 2 diverge in
both versions) are due to SciPy: since version 1.15, it uses a C
translation of MINPACK, which rounds differently than the original Fortran
code in rare cases (the first difference is in the last bit of a trial
step, i.e. in `lmpar`/`qrsolv`). With the Fortran LMDER, these fits are
identical as well. On NGC6048 this changes one Gaussian of wavelet scale 1;
the resulting residual image changes the wavelet decomposition at the 1e-14
level, which changes one Gaussian of scale 2.

## Files

| File | Purpose |
| --- | --- |
| `record_fits.py <out.pkl> <fits> [kwargs]` | Records all `lmder_fit` calls of a run with `_cbdsm` |
| `replay_fits.py <records.pkl> <mgfunction.py> [repeats]` | Compares results and speed of single fits |
| `run_comparison.sh <work dir> [rev] [repeats]` | Unit tests, end-to-end timing, and exact comparison of all FITS outputs |

The scripts need `venv` with an installed PyBDSF that still has `_cbdsm`,
as `../rmsimage-optimization/setup_overlays.sh` does. Example:

```bash
venv/bin/python -I benchmark/replace-cbdsm/record_fits.py /tmp/test.pkl test/tbdsf_process_image.in_fits \
    "{'beam': (0.015, 0.015, 0.0), 'thresh_isl': 3.0, 'thresh_pix': 5.0, 'atrous_do': True}"
venv/bin/python -I benchmark/replace-cbdsm/replay_fits.py /tmp/test.pkl bdsf/mgfunction.py
benchmark/replace-cbdsm/run_comparison.sh /tmp/cbdsm-bench master
```

# Replace the C++ extension `_cbdsm` with Python

Branch: `replace_cbdsm`, on top of `master` (7ddb350).

## What changed

- **New module `bdsf/mgfunction.py`** with `MGFunction`, `Gtype` and
  `lmder_fit`, which have the interface of those of `_cbdsm` that
  `gausfit.py` uses (the unused `dn2g_fit` and `dnsg_fit`, and the unused
  `MGFunction` methods for removing, getting and setting single Gaussians,
  were dropped). `gausfit.py` imports them
  from there; nothing else in `gausfit.py` changed.
- **Fitting:** `lmder_fit` uses MINPACK's LMDER through
  `scipy.optimize.leastsq` (the same algorithm, with the same settings, as
  the C++ code), with an analytical Jacobian.
- **Exact arithmetic:** the fits of some islands are ill-conditioned: a 1-ulp
  difference in the function or its Jacobian can lead to a different
  result. A first version (one matrix multiplication for all Gaussians,
  float64 `exp`) made single fits about 2× cheaper than the final version,
  but changed 2 of the 50 Gaussians of the test image (total flux by +67%
  and +10%), and the test failed: the PSF-variation maps differed from the
  reference outputs by 3–32%. (That version didn't cache fit results yet,
  and was about 5% slower than C++ end-to-end on the test image.) The
  Gaussians and their derivatives are therefore
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
  repeats of an earlier one (74% of the fits for the test image, 84% for
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
| test image (904) | 904 | 0.55 s | 1.87 s | 3.4 |
| NGC6048 (33097) | 33089 | 39.2 s | 47.0 s | 1.20 |

For NGC6048, the ratio is 7.7 for islands with < 50 pixels, and 0.91 for
islands with ≥ 1000 pixels.

**End-to-end** (`run_comparison.sh`, 3 interleaved runs each):

| | C++ | Python |
| --- | --- | --- |
| `test/tbdsf_process_image.py` (2 cores) | 6.59–6.81 s | 6.30–6.36 s |
| `run_pybdsf.py` NGC6048 (20 cores), total | 63.9–64.1 s | 54.4–54.9 s |
| … of which `Op_gausfit` | 12.6–12.8 s | 7.6–8.0 s |
| … of which `Op_wavelet_atrous` | 40.0–40.2 s | 35.8–36.0 s |

**Memory** (peak total PSS of the process tree, which counts pages shared
by forked workers once; interleaved runs; measured on top of `master`
4fc79a2, before the simplifications of `mgfunction`, which don't change
the peak memory of a fit):

| | C++ | Python |
| --- | --- | --- |
| `test/tbdsf_process_image.py` (5 runs) | 136–138 MiB | 136–137 MiB |
| `run_pybdsf.py` NGC6048 (6 runs) | 3730–3774 MiB (mean 3759) | 3711–3809 MiB (mean 3759) |
| … peak during `Op_gausfit` of the image | 1404 MiB | 1342 MiB |

- `gausfit.py` imports `mgfunction` at module level, so the parent process
  imports it before forking the workers. When each worker imported it, it
  copied about 2.5 MiB of shared memory pages (C++: 1.5 MiB for `_cbdsm` and
  its libraries), and the test image took 2–3 MiB more than with C++.
- The peak of NGC6048 is in the wavelet stage, in the main process, and
  doesn't depend on the fitting code; it varies by ±50 MiB between runs of
  either version.
- For the largest island of NGC6048 (2954 pixels, 10 Gaussians), a Python
  fit temporarily needs more memory than a C++ fit (4.5 MB against 1.0 MB
  on top of the memory in use before the fit, in a replay of the recorded
  fits), because SciPy copies the Jacobian that `MGFunction` returns. Python frees this memory after
  each island, while C++ keeps its caches for the largest island until the
  process ends.

- All FITS outputs (images and catalogs) are identical to those of C++:
  51 for the test image, 40 for NGC6048.

## Remaining difference

The 8 non-identical fits of NGC6048 are 2 distinct fits (one is repeated 7
times by `fit_island`), which diverge in both versions: C++ ends with NaN
parameters, Python with huge ones. Neither converges, so `gausfit` rejects
both, and the outputs are identical.

On top of `master` 4fc79a2, there were 4 non-identical fits (2 of which
diverged in both versions), which changed one Gaussian of wavelet scale 1,
and through the residual image one of scale 2. These were due to SciPy:
since version 1.15, it uses a C translation of MINPACK, which rounds
differently than the original Fortran code in rare cases (the first
difference is in the last bit of a trial step, i.e. in `lmpar`/`qrsolv`).
With the Fortran LMDER, these fits were identical as well. (This wasn't
checked again for the current 8 fits.)

## Files

| File | Purpose |
| --- | --- |
| `record_fits.py <out.pkl> <fits> [kwargs]` | Records all `lmder_fit` calls of a run, with `_cbdsm` or `mgfunction` |
| `replay_fits.py <records.pkl> <impl A> <impl B> [repeats]` | Compares results, speed and memory of single fits of two implementations (a `.py` file or a module name such as `bdsf._cbdsm`) |
| `run_comparison.sh <work dir> [rev] [repeats]` | Unit tests, end-to-end timing, and exact comparison of all FITS outputs |

`run_comparison.sh` and comparisons with C++ need `venv` with an installed
PyBDSF that still has `_cbdsm`, as `../rmsimage-optimization/setup_overlays.sh`
does. Example:

```bash
venv/bin/python -I benchmark/replace-cbdsm/record_fits.py /tmp/test.pkl test/tbdsf_process_image.in_fits \
    "{'beam': (0.015, 0.015, 0.0), 'thresh_isl': 3.0, 'thresh_pix': 5.0, 'atrous_do': True}"
venv/bin/python -I benchmark/replace-cbdsm/replay_fits.py /tmp/test.pkl bdsf._cbdsm bdsf/mgfunction.py
benchmark/replace-cbdsm/run_comparison.sh /tmp/cbdsm-bench master
```

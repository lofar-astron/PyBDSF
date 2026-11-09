# Replace the Fortran extension `_pytesselate` with Python

Branch: `replace-pytesselate`, on top of `master` (4fc79a2).

## What changed

- **New module `bdsf/tessellate.py`** with `tessellate_simple` and
  `tessellate_roundness`, which have the same interface as `pytess_simple`
  and `pytess_roundness` of `_pytesselate`. `psf_vary.py` imports it
  instead and calls these functions; nothing else in `psf_vary.py` changed.
- **Build:** `src/fortran/pytess_simple.f`, `pytess_roundness.f` and
  `constants.h` were removed, as were the f2py module in `CMakeLists.txt`
  (and `find_package(F2PY)`), the sources in `MANIFEST.in`, and `meson` from
  the build requirements (needed only by f2py). Fortran itself is still
  needed for MINPACK and PORT3 (used by `_cbdsm`). `src/fortran/gaul2gaulbin.f`
  was never built or used, and was left alone.
- **New test** `test/ttessellate.py` (added to CI and the cibuildwheel
  `test-command`).

## How it works

Each pixel belongs to the generator with the smallest (weighted) distance. The
Fortran code computes the distance to every generator for every pixel:
`n * m * ngens` distances (and twice that for `roundness`, or for fuzzy
tessellation). Here:

- The image is divided into blocks of 128 x 128 pixels. For each block, lower
  and upper bounds of the distance to each generator give the generators that
  can be the nearest one (or within the fuzzy margin) for any of its pixels.
  Most blocks have only one candidate, and so belong to that tile.
- Blocks with several candidates are divided into blocks of 8 x 8 pixels, and
  the selection is repeated (vectorized) for those, with the candidates of the
  large block.
- For the small blocks that still have several candidates, the distances are
  computed, for the candidates only, with exactly the same floating point
  operations as in the Fortran code: `sqrt((i-x)*(i-x) + (j-y)*(j-y)) / w`,
  with 1-based pixel coordinates; for `roundness`, the weight is
  `1 / (r * roundpix)`, with two roundings, as in Fortran. The nearest
  generator is selected as in the Fortran loop: the first one with a distance
  smaller than all previous ones (and than `1e90`), so ties go to the lowest
  index and NaN distances are ignored. Small blocks are evaluated in groups
  with the same number of candidates; small groups are merged by padding with
  repeats of the last candidate (which cannot change the nearest generator,
  and are excluded from the fuzzy sum).
- The bounds have a relative margin of `1e-9`, far above the rounding errors
  (a few ulp), so no generator that can be the nearest for some pixel is
  excluded. Generators for which no bounds can be given (non-finite position,
  zero, negative, infinite or NaN weight; never the case in PyBDSF) are
  candidates for every pixel; the others are still selected by their bounds.
- Small images (`n * m * ngens <= 2**20`, like the test image) are processed
  as a single large block, to avoid overhead.
- **`roundness`** first tessellates without weights, then computes the
  centroid of each tile and the mean distance of its pixels to it (the
  Fortran `tile_roundness`), and tessellates again with those weights (the
  Fortran code iterates exactly twice). The numbers of pixels and the sums of
  the coordinates are sums of integers (exact), and are collected per block
  during the first tessellation. The sum of the distances depends on the
  order of summation, so it is done in the order of the Fortran loop
  (`N.bincount` adds sequentially, and the sums are carried over from chunk to
  chunk of 16384 pixels). The result of the first tessellation is not needed
  after that, so its memory is reused for the second one.
- No array of `n * m` integers is needed (Fortran: `minind`), so the peak
  memory is the output array plus small temporary arrays.

### Deviations (undefined behaviour in Fortran)

- **`tessellate_roundness` with fuzzy tessellation (`code != 's'`)**: the Fortran
  code computes the tile properties from the fuzzy tessellation, whose values
  are sums of tile numbers, and uses them as array indices. This writes out of
  bounds: it crashed (`munmap_chunk(): invalid pointer`, `double free or
  corruption`) for both the test image and NGC6048. Here the tile properties
  are computed from the normal tessellation, which is what the code intends.
  PyBDSF does not process the result of fuzzy tessellation (see the
  documentation of `psf_tess_sc`), and the default is `'s'`.
- Pixels for which no generator has a valid distance (only possible for
  invalid input) are 0; in Fortran, they were uninitialized.
- Without generators the result is all 0 (Fortran: undefined).

## Results

Machine: 20-core i7-12700, NumPy 2.4.6. Baseline: the `_pytesselate` from the
build of `master` installed in `venv`.

**Single calls** (`bench_tess.py`, best of at least 3 calls; peak = increase of
the peak RSS during the first call). Generators: those of the recorded calls
for the test image (`test/tbdsf_process_image.in`) and for NGC6048 (the
options of `benchmark/run_pybdsf.py` with `psf_vary_do=True`), and random ones.
All results are identical to those of Fortran.

| Case | Size | Generators | Fortran | Python | Ratio | Peak MB (F / P) |
| --- | --- | --- | --- | --- | --- | --- |
| test image, simple unity 's' | 129 x 129 | 7 | 0.32 ms | 0.26 ms | 1.2 | 0 / 0 |
| test image, simple log10 's' | 129 x 129 | 7 | 0.31 ms | 0.29 ms | 1.1 | 0 / 1 |
| test image, simple sqrtlog10 's' | 129 x 129 | 7 | 0.31 ms | 0.26 ms | 1.2 | 0 / 1 |
| test image, roundness 's' | 129 x 129 | 7 | 0.95 ms | 0.83 ms | 1.1 | 0 / 1 |
| test image, simple unity 'c' | 129 x 129 | 7 | 0.82 ms | 0.39 ms | 2.1 | 0 / 1 |
| test image, simple log10 'c' | 129 x 129 | 7 | 0.84 ms | 0.44 ms | 1.9 | 0 / 1 |
| test image, simple sqrtlog10 'c' | 129 x 129 | 7 | 0.81 ms | 0.42 ms | 1.9 | 0 / 1 |
| NGC6048, simple unity 's' | 5762 x 5762 | 15 | 1.384 s | 0.146 s | 9.5 | 379 / 255 |
| NGC6048, simple log10 's' | 5762 x 5762 | 15 | 1.361 s | 0.116 s | 11.8 | 379 / 254 |
| NGC6048, simple sqrtlog10 's' | 5762 x 5762 | 15 | 1.357 s | 0.134 s | 10.1 | 379 / 255 |
| NGC6048, roundness 's' | 5762 x 5762 | 15 | 3.844 s | 0.571 s | 6.7 | 379 / 255 |
| NGC6048, simple unity 'c' | 5762 x 5762 | 15 | 3.766 s | 0.362 s | 10.4 | 379 / 255 |
| NGC6048, simple log10 'c' | 5762 x 5762 | 15 | 3.732 s | 0.252 s | 14.8 | 379 / 256 |
| NGC6048, simple sqrtlog10 'c' | 5762 x 5762 | 15 | 3.728 s | 0.323 s | 11.6 | 379 / 255 |
| random, simple log10 's' | 5762 x 5762 | 50 | 4.295 s | 0.162 s | 26.5 | 379 / 254 |
| random, simple log10 'c' | 5762 x 5762 | 50 | 12.576 s | 0.283 s | 44.5 | 379 / 255 |
| random, roundness 's' | 5762 x 5762 | 50 | 11.653 s | 0.669 s | 17.4 | 379 / 254 |
| random, simple, w[0] = 0, 's' | 5762 x 5762 | 50 | 4.361 s | 0.718 s | 6.1 | 380 / 255 |
| random, simple log10 's' | 5762 x 5762 | 200 | 17.209 s | 0.264 s | 65.2 | 379 / 254 |
| random, simple log10 'c' | 5762 x 5762 | 200 | 50.713 s | 0.412 s | 123 | 379 / 255 |
| random, roundness 's' | 5762 x 5762 | 200 | 45.713 s | 0.990 s | 46.2 | 379 / 255 |
| random, simple, w[0] = 0, 's' | 5762 x 5762 | 200 | 17.439 s | 0.757 s | 23.0 | 379 / 255 |
| random, simple log10 's' | 2000 x 2000 | 1000 | 10.233 s | 0.089 s | 115 | 45 / 32 |
| random, simple log10 'c' | 2000 x 2000 | 1000 | 30.395 s | 0.140 s | 217 | 45 / 32 |
| random, roundness 's' | 2000 x 2000 | 1000 | 27.020 s | 0.302 s | 89.5 | 46 / 33 |
| random, simple, w[0] = 0, 's' | 2000 x 2000 | 1000 | 10.270 s | 0.138 s | 74.3 | 45 / 32 |

The NGC6048 rows are from a separate run of `bench_tess.py` with only those
cases: during the full run, other processes on the machine slowed down the
Python version (by up to 3x); the other rows are from the full run. The
Fortran version of `roundness` with code 'c' crashed in both cases (see
above).

**End-to-end** (`run_comparison.sh`):

| | Fortran | Python |
| --- | --- | --- |
| `test/tbdsf_process_image.py` | passed | passed |
| NGC6048 (`run_psf_vary.py`), simple tessellation | 1.21 s | 0.12 s |
| NGC6048 (`run_psf_vary.py`), roundness tessellation | 3.40 s | 0.50 s |

- Test image: all 48 FITS outputs (images and catalogs) are identical.
- NGC6048: all 44 FITS outputs (including the `volrank` image and the PSF
  maps) are identical, for both tessellation methods.
- The total run time of NGC6048 (about 70-90 s) varied more between runs than
  the time saved.

## Notes

- **Small images and the heap:** in a process with a small heap, glibc
  returns the memory freed after each call to the system, and the next call
  takes page faults, which made the Python version up to 30% slower on the
  test image (1.1-1.2 ms instead of 0.85 ms for `roundness`). This does not
  happen in a PyBDSF run (or with `MALLOC_TRIM_THRESHOLD_` and
  `MALLOC_MMAP_THRESHOLD_` set high), so `bench_tess.py` imports `bdsf` in the
  worker processes of both versions.
- **Clean build:** a wheel built from a clean copy of this branch
  (`pip wheel .`, with build isolation) contains `tessellate.py` and no
  `_pytesselate`, and `tmulti_proc.py` and `ttessellate.py` pass with it.
  `tbdsf_process_image.py` crashes in the Gaussian fitting (`_cbdsm`) with
  such a wheel, but in exactly the same way with a wheel of `master` built in
  the same way, so this is unrelated.
- The reference implementation in `test/ttessellate.py` (a direct
  transcription of the Fortran loops) was checked against `_pytesselate`:
  identical in all 216 random and degenerate cases (ties, duplicate
  generators, empty tiles, tiles of one pixel), and in all 22 cases with
  invalid weights or positions (where the Fortran result is defined).

## Files

| File | Purpose |
| --- | --- |
| `record_tess.py <out.pkl> <image> [kwargs]` | Records all `_pytesselate` calls of a PyBDSF run |
| `bench_tess.py <tessellate.py> [records.pkl ...] [--synthetic] [--repeat N]` | Compares results, speed and peak memory of single calls |
| `run_psf_vary.py <image> [psf_itess_method]` | Runs NGC6048 with PSF variation, and reports the time of the tessellation |
| `compare_fits.py <dir1> <dir2>` | Checks that all FITS outputs of two runs are identical |
| `run_comparison.sh <work dir> [rev]` | Unit test, test image and NGC6048 with both versions, with exact comparison of all FITS outputs |

The scripts need `venv` with an installed PyBDSF that still has
`_pytesselate`, as `../rmsimage-optimization/setup_overlays.sh` does. Example:

```bash
venv/bin/python -I benchmark/replace-pytesselate/record_tess.py /tmp/ngc.pkl benchmark/NGC6048_9-MFS-image.fits \
    "{'thresh_isl': 4.0, 'thresh_pix': 5.0, 'rms_box': (150, 15), 'adaptive_rms_box': True, \
      'adaptive_thresh': 150, 'rms_box_bright': (60, 15), 'frequency': 143.65e6}"
venv/bin/python -I benchmark/replace-pytesselate/bench_tess.py bdsf/tessellate.py /tmp/ngc.pkl --synthetic
benchmark/replace-pytesselate/run_comparison.sh /tmp/tess-bench master
```

(`record_tess.py` writes the outputs of the run next to the image, so it is
best run on a copy or a symbolic link in a separate directory.)

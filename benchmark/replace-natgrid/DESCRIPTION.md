# Replace the C extension `natgridmodule` with Python

Branch: `replace-natgrid-extension`, on top of `53d19e8`.

## What changed

- **New module `bdsf/natgrid.py`** with one function,
  `natgrid(x, y, z, xo, yo)`: natural neighbor (Sibson) interpolation of
  scattered data to a rectangular grid, as `natgridd` of NCAR's `natgrid`
  (Dave Watson's `nngridr`) does with its default parameters, which is all
  PyBDSF uses. `z` may also be 2-D, to interpolate several sets of values
  with the same points at once.
- **`psf_vary.py`**: `interp_prop` and `interp_shapcoefs` call `natgrid`
  instead of `nat.Natgrid(...).rgrd(...)`. `interp_shapcoefs` now interpolates
  all shapelet coefficients in one call.
- **Removed**: the `natgrid/` sources and the `bdsf.nat` package (the
  object-oriented wrapper).
- **Build**: there is no compiled code left (`src/fortran/gaul2gaulbin.f` was
  never built), so PyBDSF is now built as a pure Python package with
  setuptools: `CMakeLists.txt`, `setup.py`, the scikit-build, CMake, Ninja and
  NumPy build requirements, the cibuildwheel settings and `cibuildwheel/` were
  removed (an empty CMake project cannot be built with scikit-build: it has no
  `install` target), as were the C and Fortran classifiers and the compiler
  notes in `README.rst`. The release workflow builds one `py3-none-any`
  wheel and the sdist with `build`, and runs the tests with the wheel
  (cibuildwheel refuses to build pure Python wheels).
- **New test** `test/tnatgrid.py`, with reference outputs of the C extension
  in `test/natgrid_reference.npz` (made with `make_reference.py`), added to
  CI and the cibuildwheel `test-command`.
- **Updated test reference outputs**: the `psf_vary_*` images in
  `test/reference_outputs/misc` (see the bug below).

## Bug in the old code: scrambled tile coordinates

`psf_vary.py` passes `x, y = N.transpose(psfcoords)` to `nat.Natgrid`. These
are strided views (every other element of the `(n, 2)` array), and
`natgridmodule` uses `PyArray_DATA` without checking for contiguity, so it
read `x = [x0, y0, x1, y1, ...]` and `y = [y0, x1, y1, x2, ...]` (the first n
values from the start of each view). The PSF variation maps
(`psf_vary_maj/min/pa/ratio` and the aperture ratio), and therefore the PSF
used for deconvolved sizes, were interpolated from tile positions that have
nothing to do with the actual tiles. (The comment "else natgrid cant deal
with noncontiguous memory" next to the shapelet coefficients shows the
problem was known for the values, but not for the coordinates.) This was
verified by calling `natgridmodule` with these scrambled coordinates: the
result is bitwise identical to the output in a real run.

`bdsf.natgrid` reads its input with `numpy.asarray`, so it uses the actual
coordinates. Its results can therefore not be compared directly with the old
outputs; `run_comparison.sh` compares them with those of the old code with
contiguous coordinates (`fix`).

## How it works

Like `natgrid` (`Initialized` and `ReadDatad`):

- The overlap region (the output grid extended to include all data points) is
  determined, with the limits rounded to single precision as in `armind` and
  `armaxd`; only data points inside it are used.
- If the extents of the region differ by more than a factor 2, an axis is
  rescaled (automatic scaling).
- The three vertices of a large enclosing triangle (60 times the extent) are
  added, with the values of the least-squares plane through the data. These
  determine the interpolant outside the convex hull of the data
  (extrapolation) and close to it.
- The data coordinates are perturbed by up to `1e-5 / 2` of the extent, with
  the `rand()` sequence after `srand(367)`. The glibc sequence is reproduced,
  so the results are the same as those of the C code on Linux (and now also
  the same on macOS).
- Duplicate data points raise a `ValueError` (natgrid error 2).

The interpolation itself is different. `natgrid` loops, for every grid
point, over all Delaunay triangles, and for each triangle whose circumcircle
contains the grid point, computes the circumcenters of the grid point and each
edge, and adds the signed areas of the triangles they form with the
circumcenter of the triangle to the weights of the vertices (Watson's method).
These areas have a closed form: with `A` the area of the triangle, `l0, l1, l2`
the barycentric coordinates of the grid point q and `D = R^2 - |q - c|^2`
(positive inside the circumcircle with center c and radius R), the area added
to vertex k is `D^2 l_k / (16 A l0 l1 l2)`. So a triangle contributes
`W = D^2 / (16 A l0 l1 l2)` to the total weight, and `W * f(q)`, with f the
linear interpolant over the triangle, to the weighted sum. The interpolated
value is the ratio of the two sums.

- The Delaunay triangulation (of the data and the enclosing triangle) is
  computed with `scipy.spatial.Delaunay`.
- The grid is processed in blocks of 32 rows. For each block and triangle,
  `W` and `W * f` are evaluated (with in-place NumPy operations on
  preallocated buffers) only on the rows and columns that the circumcircle
  covers, and are zero outside it.
- Close to an edge, the contributions of the two triangles that share the
  edge are large and nearly cancel. To make their rounding errors cancel as
  well, the barycentric coordinate of the opposite vertex is computed from the
  line function of the edge in exactly the same way in both triangles (same
  vertex order and origin). `D` and `f` are computed relative to a data
  vertex of the triangle, not a vertex of the enclosing triangle. Without
  this, errors up to `2e-2` occurred along such lines in the NGC6048 run.
- At grid points exactly on an edge, the contributions are infinite; such
  points (so far never seen with real data, since the coordinates are
  perturbed) get the value at a point shifted by `1e-9` of the extent.
- Memory: the output array plus buffers for one block (a few MB). The C code
  also allocated the full output array.

### Accuracy

For the NGC6048 interpolations, the maximum difference with the C code is
`1.8e-7` (on values of about 2.5); only 450 of 33 million pixels differ by
more than `1e-10`, all along lines through an edge. Compared with exact
(rational) evaluation of the same interpolant at the worst pixels, the error
of the C code is up to `1.6e-7` and that of the Python code up to `1.2e-8`;
at 200 random pixels the maximum errors are `1.4e-12` (C) and `1.0e-13`
(Python).

### Deviations

- Decreasing output coordinates are supported (the result is reversed); the
  C extension crashed with them (segmentation fault).
- The output is always `float64` (the C code used single precision for
  `float32` input; PyBDSF always passes `float64`).
- Only the default control parameters are supported (no gradients, aspects or
  slopes, single-point mode, wrapping, or settable parameters); PyBDSF used
  none of these.

## Results

Machine: 20 cores, numpy 2.5.5, scipy 1.18.1. `run_comparison.sh <dir> HEAD 3`.

### Single interpolation (NGC6048, 13 tiles, 5762 x 5762 grid)

One of the four `interp_prop` calls, recorded with `record_nat.py` and timed
with `bench_call.py`, each in a fresh process (the other calls take the same
time):

| version | CPU time | peak RSS increase |
|---------|---------:|------------------:|
| C       | 3.41 s   | 253 MB (= output) |
| Python  | 1.10 s   | 263 MB (output + 10 MB of buffers and scipy) |

(Directly after `import numpy`, OpenBLAS threads spin for about half a second;
measurements must not start in that window, or set `OPENBLAS_NUM_THREADS=1`.)

### `run_pybdsf.py` (NGC6048, `psf_vary_do=True`)

| version | CPU time (user + sys) | peak RSS | `Op_psf_vary` (wall) | total (wall) |
|---------|----------------------:|---------:|---------------------:|-------------:|
| orig    | 279.7 / 278.3 / 277.6 s | 5709 / 5698 / 5664 MB | 7.97 / 6.98 / 6.78 s | 68.3 / 65.2 / 65.4 s |
| fix     | 281.9 / 275.2 / 277.5 s | 5686 / 5696 / 5679 MB | 8.33 / 6.97 / 6.96 s | 69.0 / 65.0 / 65.0 s |
| new     | 275.4 / 269.6 / 267.3 s | 5719 / 5671 / 5719 MB | 6.68 / 4.73 / 4.70 s | 69.1 / 63.1 / 63.0 s |

`Op_psf_vary` is about 2.2 s faster (the four interpolations run in parallel,
so this is the time of about one of them), and the CPU time is about 9 s
lower. The peak memory is the same (it is reached elsewhere, and the
interpolation needs the same memory as the C code).

### Outputs

- `run_pybdsf.py`, new vs fix: 40 of 46 FITS files identical; the
  `psf_vary_*` images differ by at most `2e-7` (relative; these are single
  precision images, so this is about one ulp), and the catalogs by at most
  `3e-12`.
- `run_pybdsf.py`, new vs orig: the `psf_vary_*` images differ by up to 130%,
  and the deconvolved sizes in the catalogs differ, because of the bug fix.
- `test/tbdsf_process_image.py`, new vs fix: 45 of 51 FITS files identical;
  the `psf_vary_*` images differ by at most `3.3e-12`, the catalogs by at most
  `5e-14`. With the old reference outputs, both fix and new fail (on the
  `psf_vary_min/pa/ratio` images); the updated references are those of new.
- PSF shapelet interpolation (`psf_use_shap=True`, test image): `psf_images`
  identical for new and fix.
- All tests (`tgausfit_shape`, `tmgfunction`, `tmulti_proc`, `tnatgrid`,
  `ttessellate`, `tbdsf_process_image` with default and CI tolerance) pass
  with the wheel built by `python -m build`, installed in a fresh venv.

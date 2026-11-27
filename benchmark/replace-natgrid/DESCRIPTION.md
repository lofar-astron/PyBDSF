# Replace the C extension `natgridmodule` with Python

The natgrid C extension (NCAR ngmath, used by `psf_vary.py`) is replaced by
`bdsf/natgrid.py`. How the interpolation works is described in its module
docstring.

## Files

- `run_comparison.sh <work dir> <baseline revision> [repeats]`: runs the unit
  test, `test/tbdsf_process_image.py` and `benchmark/run_pybdsf.py` with
  three versions, and compares their FITS outputs:
  - `orig`: the baseline, with the C extension;
  - `fix`: the baseline, with contiguous coordinates passed to the C
    extension (see the bug below); the new code is compared with this one;
  - `new`: the working tree, without the C extension.
- `record_nat.py`: runs `run_pybdsf.py` and saves the inputs and outputs of
  every natgrid call.
- `bench_call.py`: CPU time and peak memory of one recorded call, with the C
  extension or with `bdsf/natgrid.py`.
- `make_reference.py`: creates `test/natgrid_reference.npz` (outputs of the
  C extension, used by `test/tnatgrid.py`).

The scripts that use the C extension need a PyBDSF installation that still
has `bdsf.nat`.

## Bug in the old code: scrambled tile coordinates

`psf_vary.py` passed `x, y = N.transpose(psfcoords)`, which are strided views,
to `natgridmodule`, which ignores strides. It read `x = [x0, y0, x1, y1, ...]`
and `y = [y0, x1, y1, x2, ...]`, so the PSF variation maps were interpolated
from wrong tile positions. Calling the C module with these scrambled
coordinates reproduces the output of a real run bit for bit.
`bdsf/natgrid.py` uses the actual coordinates, so its results are compared
with `fix`, and the `psf_vary_*` reference images of
`test/tbdsf_process_image.py` were regenerated.

## Differences from the C code

- The results agree up to rounding errors. Near lines through edges of the
  triangulation, where the contributions of two triangles nearly cancel, the
  Python code is more accurate: compared with exact (rational) evaluation at
  the pixels where the two differ most, the errors are up to `1.6e-7` (C)
  and `1.2e-8` (Python).
- Only natgrid's default parameters are supported, and the output
  coordinates must be increasing; this is all PyBDSF uses. (The C extension
  crashed with decreasing output coordinates.)
- The output is always `float64`.

## Results

Machine with 20 cores, numpy 2.5.5, scipy 1.18.1;
`run_comparison.sh <dir> 53d19e8 3` and `bench_call.py`.

One interpolation (NGC6048: 13 tiles, 5762 x 5762 grid), in a fresh process:

| version | CPU time | peak RSS increase |
|---------|---------:|------------------:|
| C       | 3.37 s   | 253 MB (the output) |
| Python  | 1.10 s   | 263 MB |

(Right after `import numpy`, OpenBLAS threads spin for about half a second;
`bench_call.py` waits 2 s before it starts timing.)

`run_pybdsf.py` (NGC6048, `psf_vary_do=True`), 3 runs each:

| version | CPU time (user + sys) | peak RSS | `Op_psf_vary` (wall) | total (wall) |
|---------|----------------------:|---------:|---------------------:|-------------:|
| orig | 286.4 / 279.2 / 278.4 s | 5669 / 5697 / 5672 MB | 9.13 / 6.76 / 6.84 s | 71.8 / 65.3 / 65.0 s |
| fix  | 285.2 / 278.5 / 278.9 s | 5663 / 5735 / 5670 MB | 9.70 / 7.01 / 6.96 s | 72.3 / 65.0 / 65.8 s |
| new  | 279.1 / 268.7 / 268.1 s | 5673 / 5680 / 5728 MB | 7.10 / 4.63 / 4.61 s | 69.4 / 62.5 / 63.4 s |

The peak memory is reached elsewhere in the pipeline, and is unchanged.

Outputs, new vs fix:

- `run_pybdsf.py`: 40 of 46 FITS files identical; the `psf_vary_*` images
  differ by at most `2e-7` (relative; about one float32 step), the catalogs
  by at most `7e-12`.
- `test/tbdsf_process_image.py`: 45 of 51 FITS files identical; the others
  differ by at most `3.3e-12`.
- PSF shapelet interpolation (`psf_use_shap=True`, test image): identical
  `psf_images`.

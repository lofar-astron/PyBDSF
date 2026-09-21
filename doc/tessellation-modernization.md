# Tessellation modernization work log

## Check-script reporting update

The equivalence script now reports reference setup, progress through deterministic,
randomized and large-image cases, invalid-input and intermediate-roundness checks,
per-mode comparison totals, and elapsed time. Use `--verbose` (or `-v`) to see
every comparison's name, image dimensions, generator count, mode, epsilon and
result. Failed comparisons always print their context before the traceback.
The numerical checks and generated inputs are unchanged.

## Objective

Replace the two tessellation routines in `src/fortran` with Python source using
the existing NumPy/SciPy dependencies, provided that defined legacy outputs are
identical and representative performance is at least as good. No JIT compiler
or new native extension is planned. Keep the original sources as an independent
test and benchmark reference.

## Inspection (2026-09-21)

- `bdsf/psf_vary.py` calls only `pytess_simple` and `pytess_roundness`.
- Coordinates and output labels are one-based; output arrays are float64.
- Exact ties go to the first generator. Simple distances are evaluated as
  `sqrt(dx*dx + dy*dy) / weight`, not squared distances.
- Simple fuzzy output is an integer-valued sum of pairs of generator labels,
  not a nearest-generator label map.
- Roundness runs exactly two assignment passes. The first is unweighted; the
  second uses distances to the first-pass centroids and mean radii.
- Roundness reductions traverse columns first (Fortran memory order). Preserving
  this order is necessary for exact floating-point equivalence.
- Fuzzy roundness is unsafe in the reference: fuzzy values can be zero or
  greater than the generator count, but are used as array indices. There is no
  meaningful general output-equivalence contract for this mode.
- Other Fortran libraries (`minpack`, `port3`) still serve the C++ fitting
  extension. Replacing tessellation alone will not remove the Fortran compiler
  requirement for the whole project.

## Approach under evaluation

Use vectorized NumPy arithmetic and conservative spatial bounds to discard
generators that cannot win a rectangular block. Evaluate retained candidates in
their original order and preserve legacy arithmetic at pixel boundaries.
Measure against a separately compiled, optimized copy of the original routines.

Implementation, correctness results, performance measurements, and commands will
be recorded below as the work proceeds.

## Implementation and experiments

- Added `bdsf/_tessellation.py`, using NumPy only. There are no new dependencies,
  compiled modules, JIT compilation, approximate nearest-neighbor searches, or
  parallel worker pools.
- The initial fixed-block implementation matched the reference but was slower
  for small generator counts. Batching 16-by-16 blocks, writing directly into
  the Fortran-order output, and retaining only possible candidates substantially
  reduced its overhead.
- Bounds include a conservative relative margin. Actual decisions use the
  original square root and division expressions, and strict comparisons in
  generator order. The two divisions in roundness are intentionally retained.
- Centroids use exact integer sums of constant-label runs within each column.
  Above the float64 exact-integer range, the original accumulation order is
  used instead. Mean radii always accumulate in original column-first order.
- Single-generator cases return the constant result directly. Very small maps
  use a bounded vectorized distance calculation instead of spatial pruning.
- `psf_vary.py` now imports `_tessellation`. This distinct name prevents an old
  `_pytesselate` shared library from silently overriding the Python code.
- Removed only the tessellation F2PY target and F2PY discovery from CMake.
  The original Fortran files remain unchanged for independent verification.
- Performance scope confirmed with the user: prioritize images of 512 by 512
  and larger. Tiny images can remain slower because of Python dispatch and
  validation overhead; no universal speed guarantee is claimed.

## Reproducible checks

From the repository root, using an environment with NumPy:

```sh
python test/check_tessellation.py --random-cases 1000 --large-cases 8
python test/benchmark_tessellation.py --repeat 7 --min-speedup 1 \
    --output doc/tessellation-benchmark.json
```

The first run builds an independent Fortran reference with `gfortran -O3
-fno-fast-math` through F2PY. This requires a C compiler and gfortran; NumPy's
Meson backend also requires Meson/Ninja. Compiled references are cached under
the system temporary directory using the source, compiler, Python, NumPy and
flag versions as the cache key. See the cached `build.log` on compilation errors.
PyBDSF itself does not need to be built or installed to run these scripts.

The equality script checks every element, dtype, shape and output storage order,
with **no numerical tolerance**. It covers weighted and unweighted simple
tessellation, fuzzy simple tessellation, roundness, exact and near ties, duplicate
generators, empty and singleton tiles, off-image generators, one-dimensional
images, block boundaries, weight extremes, randomized inputs and larger images.
It separately compares roundness centroids and inverse mean radii with the
original Fortran helper, including empty-tile NaNs.

Final result: **11,236 exact output comparisons passed**, along with 10
invalid-input checks and exact intermediate roundness reductions (seed 20260921).
This is evidence over those inputs, not a proof of equivalence for all float64
inputs or compiler configurations.

The benchmark checks equality first, warms both implementations, alternates
measurement order, includes output allocation, and reports medians and every
sample. Compilation, imports, input generation and checks are excluded from
timing. `--sizes`, `--generators`, `--repeat` and `--seed` are configurable.
`--min-speedup 1` gives a nonzero exit if any case is slower than Fortran.

## Measured performance

Host: Intel Core i7-12700, Linux x86_64, Python 3.12, NumPy 1.26.4,
gfortran 13.3. Seven repeats per case, seed 20260921. No extra worker threads
are started by the implementation.

All 32 combinations of sizes 512 and 1024, generator counts 2, 8, 64 and 256,
and unity/weighted/fuzzy/roundness modes passed exact equality and the
`--min-speedup 1` gate. Observed median speedups range from **1.05x to 41.01x**.
Raw timings and environment metadata are in `tessellation-benchmark.json`.

A separate 2048 x 2048 run (2, 8 and 64 generators, all four modes, three
repeats) also passed every equality and performance check: **1.33x to 30.03x**.
Raw samples are in `tessellation-benchmark-2048.json`. Reproduce with:

```sh
python test/benchmark_tessellation.py --sizes 2048 --generators 2 8 64 \
    --repeat 3 --min-speedup 1 --output doc/tessellation-benchmark-2048.json
```

Both final reports record the implementation's SHA-256 so that timings can be
associated with the exact source that was measured.

| Image | Generators | Mode | Fortran (ms) | Python (ms) | Speedup |
| --- | ---: | --- | ---: | ---: | ---: |
| 512 x 512 | 2 | Roundness | 6.634 | 5.679 | 1.17x |
| 512 x 512 | 8 | Weighted | 5.497 | 1.722 | 3.19x |
| 1024 x 1024 | 2 | Roundness | 27.035 | 25.636 | 1.05x |
| 1024 x 1024 | 64 | Weighted | 173.638 | 14.284 | 12.16x |
| 1024 x 1024 | 256 | Fuzzy simple | 2056.478 | 50.140 | 41.01x |
| 1024 x 1024 | 256 | Roundness | 1869.819 | 112.489 | 16.62x |

These are measured workload results, not a guarantee on other machines or
generator layouts. In particular, the narrowest margin (two-generator
roundness at 1024 x 1024) should be rechecked on deployment hardware. Tiny-image
exploratory runs were slower despite exact outputs; they are outside the
user-confirmed performance target.

## Build verification

Built a complete PyBDSF wheel successfully from the modified source with
`pip wheel --no-deps --no-build-isolation`. The remaining C++/Fortran fitting
and natgrid extensions still build; the tessellation F2PY target is absent.
Temporary build dependencies and products were kept outside the user's checkout.
Installed the wheel in an isolated environment and verified that both
`Op_psf_vary.tess_simple` and `Op_psf_vary.tess_roundness` run through the packaged
Python module. Inspected the wheel contents: `_tessellation.py` is present and
there is no `_pytesselate` shared library.
The installed wheel also passed the repository's existing
`test/tbdsf_process_image.py` end-to-end test, including the PSF-variation stage
and comparison against reference outputs. Its multiprocessing workers required
execution outside the filesystem/socket sandbox.
After adding the final singleton-at-block-boundary regression case, the exact
comparison suite, wheel build, installed end-to-end test, and both benchmark
matrices were rerun successfully. Ruff and `git diff --check` also passed.

## Defined behavior and limitations

- Supported numerical inputs have positive image dimensions, nonempty finite
  coordinate arrays, matching SNR/weight array lengths, finite `eps`, positive
  finite simple weights, and at least one usable generator distance below the
  legacy `1e90` sentinel at every pixel. SNR values are unused, as in Fortran.
- Both hard and fuzzy simple outputs preserve the original encoding exactly.
  Roundness supports hard labels. Fuzzy roundness raises a clear `ValueError`
  instead of reproducing undefined memory accesses. The existing option help
  already advises against using fuzzy results downstream.
  A separate bounds-checked gfortran probe confirmed this: a 2 x 2 image with
  one generator and code `'c'` terminates at `pytess_roundness.f:152` because
  index zero is below the lower bound of the centroid array. This is a legacy
  defect, not a valid output that the replacement can reproduce.
- Empty or zero-radius roundness tiles cannot win the second pass, just as their
  NaN/infinite reference distances cannot satisfy a strict minimum comparison.
  If no usable tile remains (for example a one-pixel roundness map), Python
  raises `ValueError` instead of returning an uninitialized Fortran label.
- Validation intentionally rejects invalid inputs rather than emulating
  uninitialized memory or out-of-bounds accesses.
- Working memory is linear in the image size plus bounded batches of candidate
  bounds/distances; there is no image-sized cube for every generator. Roundness
  still requires several full-image arrays for exact reductions.
- The original private F2PY helper entry points are not a supported application
  API. The two public-to-the-caller routine signatures (including optional
  `ngens`) are preserved. The Fortran helper routines remain in the test oracle.

## Comparison with the alternative Python port

The implementation in branch `python-tessellation-opencode-big-pickle`, inspected
in `/home/marcel/code/PyBDSF.opencode/bdsf/tess.py`, reproduced the reported
`duplicate-generators/roundness` failure: **390 of 391 pixels differed** from
Fortran. This investigation did not modify either implementation.

The case uses a 17 x 23 image and generators at `(4, 6)`, `(4, 6)` and
`(15, 19)`. Ties go to the first generator, so the first assignment gives:

| Generator | First-pass pixels |
| --- | ---: |
| 1 | 212 |
| 2 | 0 |
| 3 | 179 |

In the alternative port's `_tile_roundness`, masked divisions leave the empty
tile's centroid and inverse mean radius (`roundpix`) at zero. The second-pass
score multiplies by `roundpix**2`, making that generator's score zero throughout
the image. It consequently captures 390 pixels; generator 1 retains one pixel
through tie-breaking. The Fortran result instead assigns 211 pixels to generator
1 and 180 to generator 3. The Codex implementation matched that result exactly.

Fortran computes NaN quantities for the empty tile. Its strict
`distance < minimum` comparison never selects that candidate. The Codex port
explicitly excludes unusable candidates. A correction to the alternative port
should exclude empty and zero-radius tiles from the second-pass minimum search;
substituting zero for an undefined inverse radius changes the assignment rule.

Two additional compatibility concerns were identified by inspection, separately
from this reproduced failure: the alternative port substitutes squared-distance
expressions and accumulates radii in C order rather than the original Fortran
order. These can change near-tie results through floating-point rounding. Fixing
the empty-tile bug alone therefore does not establish exact equivalence.

## Why NumPy can outperform the original Fortran

The main advantage is a reduction in the amount of work. The original Fortran
compares every pixel with every generator. The Python implementation evaluates
conservative distance bounds for 16 x 16 blocks, discards candidates that cannot
win anywhere in a block, and computes pixel-level distances only for the
remaining candidates. A block with a single candidate can be filled directly.
NumPy executes its array arithmetic in compiled native code; Python coordinates
these operations rather than executing each pixel calculation itself.

An instrumented weighted example used a 1024 x 1024 image, 256 generators and
`numpy.random.default_rng(20260921)`. Generator coordinates were drawn with
`rng.uniform(1, 1024, (2, 256))`, followed by weights from
`rng.uniform(0.5, 2.0, 256)`. Counting candidates with the implementation's bounds
and conservative margin gave:

| Work | Count |
| --- | ---: |
| Original pixel-generator distance evaluations | 268,435,456 |
| Block-generator bound pairs | 1,048,576 |
| Remaining pixel-generator distance evaluations | 1,748,736 |
| Blocks filled without pixel-level distance calculations | 1,580 of 4,096 |
| Mean retained candidates per block | 2.0535 |

Bounds cost arithmetic too, so these counts are not themselves a speedup ratio.
They explain why a measured speedup around 20x is plausible even after Python
overhead, bounds calculations and memory traffic. This count is a separate
reproducible illustration, not the exact generator sample used in every
benchmark row. The benchmark reference is optimized Fortran compiled with
`-O3 -fno-fast-math`.

Fortran using the same pruning strategy could plausibly match or outperform the
NumPy implementation. The observed gain is primarily algorithmic, not evidence
that Python arithmetic is faster than compiled Fortran. The benefit depends on
generator count and layout; with few generators or ineffective pruning, the
advantage shrinks. Worst-case pixel-distance work can still approach the
original all-pixels/all-generators search.

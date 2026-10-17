# Conversation log: release notes and rmsimage performance

Date: 2026-10-17. Participants: Marcel Loose (user) and Claude Code (Claude Opus 5.5).
Repository: PyBDSF, branch `master` at `c6efae8`.

This is a record of the conversation. Tool calls and outputs are summarised;
the key numbers are given exactly as measured.

---

## 1. Release notes for the last five commits

**User:** Can you write up a summary of the changes in the last five commits
for the release notes? The team will review it in Monday's meeting.

**Claude:** Created a shared Claude Doc, "PyBDSF Release Notes: Recent
Changes"
(https://claude.ai/code/artifact/d993f7a4-940f-4076-97ab-dc8e1fae8135),
built from `git log -5 --stat` and the diffs. Contents:

- **Overview**: a table of the five commits:

  | Commit | PR | Change | Type |
  | --- | --- | --- | --- |
  | c6efae8 | #398 | Fix WCS position error scaling for non-square pixels | Bug fix |
  | 497885a | #471 | Skip `isl_tosplit` when `split_isl` is disabled | Performance |
  | ac15748 | #467 | Warn about large islands before fitting | Usability |
  | 0cbb0b2 | #399 | Fix weighted sum of squares for the PSF position angle | Bug fix |
  | 73d70d4 | #397 | Fix underestimation of `deconv_size_skyE` | Bug fix |

- **Bug fixes**
  - #397: `deconv_size_skyE` was copied from the observed size errors. It is
    now propagated as δθ_dec = δθ_obs · θ_obs / θ_dec, for each axis
    separately. Point-like sources keep the observed error, and the PA error
    is unchanged.
  - #399: `stdpa` in `psf_vary.py` used `bpa` instead of `bpa²` in the
    weighted sum of squares.
  - #398: source `E_RA`/`E_DEC` were scaled by the geometric mean of the two
    pixel sizes; each axis now uses its own `acdelt`. (Checked: the old
    `cdeltsq` is `acdelt[0]*acdelt[1]`, so square pixels give identical
    results.)
- **Improvements**
  - #467: a warning before fitting when islands exceed 2000 active pixels
    (hard-coded threshold).
  - #471: no `isl_tosplit` call when `split_isl=False`; same output, faster.
- **Impact and discussion points**: catalogue error columns can change for
  users who upgrade. There's a checklist for the meeting: a behaviour-change
  note, making the 2000-pixel threshold configurable, the version number, and
  whether Gaussian-level position errors need the same non-square-pixel fix.

The user interrupted before the optional review comment was added to the
doc. All sections were filled in.

---

## 2. Performance of `rmsimage.py`

**User:** We run into some serious performance bottle-necks. I've traced it
down to the modules rmsimage.py. Please check the python profile file
benchmark/filter_skymodel.prof. Can you improve the performance?

### 2.1 Profile analysis (`analyze_profile.py`)

The profile covers 1007.5 s in total. The relevant entries:

| Function | cumtime | tottime | Note |
| --- | --- | --- | --- |
| `Op_rmsimage.__call__` | 703 s | 51.7 s | 5 calls |
| `rms_mean_map` | 591 s | 39.9 s | 14 calls, 42 `parallel_map` calls |
| `multi_proc.parallel_map` | 568 s | | 46 calls |
| `posix.waitpid` | | 377.6 s | parent waiting for worker processes |
| `posix.read` / pickle loads | | 95 s / 10 s | results via manager queue |
| `numpy.array` in `parallel_map` | | 31.1 s | task list converted to an object array |
| `posix.fork` | | 34.5 s | 782 forks |
| `list.append` in `rms_mean_map` | | 4.4 s | 61.7M calls |
| `ndimage.label` + `find_objects` in `__call__` | 31 s | | inherent |

Observations:

- The profile covers only the parent process. Most of the time is spent
  waiting for workers (`waitpid`), whose own work isn't visible.
- Each `rms_mean_map` call built ~2.2M boxes (61.7M appends / 14 calls / 2
  lists), so the production images are much larger than the 5762² test image.
- `rms_mean_map` called `parallel_map` 3 times (interior, two borders) with
  `zip(itertools.repeat(arr), ...)`, i.e. one tuple per box that included
  the whole image. `parallel_map` then turned that list into a numpy object
  array, forked ~17 processes per call, and each worker sent back a list of
  `(mean, rms)` tuples of numpy scalars through a `multiprocessing.Manager`
  queue.

### 2.2 First attempt: vectorized bstat (rejected)

Claude first wrote `bstat_batch`: it sorts many boxes at once, keeps the
clipping state as index bounds plus cumulative sums, and uses a row-wise
binary search. Benchmarks on the 5762² image (`rms_mean_map`, box (150, 15)):

- original: 13.4 s; `bstat_batch`: 19.2 s (**slower**).

Measuring per box showed why: the original `bstat` takes 1.20 ms per
150×150 box, of which **1.07 ms is the `sort`**. The clipping iterations are
only ~10% of the cost, and the batch version added float64 cumulative sums
and cache pressure. The approach was reverted.

### 2.3 Lower multiprocessing overhead in `rms_mean_map`

- Box indices and map coordinates are built with numpy (`meshgrid`,
  `column_stack`) instead of Python loops.
- Interior boxes, borders and the corner box go to **one** `parallel_map`
  call instead of three, with the same order of operations: the "no unmasked
  regions" check still sees only the interior boxes.
- The new method `process_mean_rms_boxes()` passes only chunk numbers to
  `parallel_map`. Workers inherit the image through fork, loop over their
  boxes with the existing `for_masked_mp()`, and return a numpy array per
  chunk.
- There are `100 * nproc()` chunks. A first version with `4 * nproc()`
  chunks was slower (67 s vs 56 s on the tiled image) because
  `parallel_map` hands each process a whole number of chunks: 80 chunks over
  19 processes gives some 5 and others 4.
- `process_mean_rms_maps()` and `for_masked()` became unused and were
  removed, as was `import itertools`.

With the original `bstat`, the maps were **bit-identical** to before. The
speed gain was small on this machine: 13.4 → 13.7 s (noise) on the 5762²
image, and 60.7 → 58.0 s on the 2×2-tiled image. It matters more for the
much larger production images in the profile.

### 2.4 `functions.bstat` without a full sort

`bstat` only needs, at each iteration, the median of the remaining pixels,
the pixels in `[med - κσ, med + κσ)`, and their mean and sum of squares.
New implementation:

- The median comes from `ndarray.partition(npix//2)` (O(n)) instead of a
  full sort.
- Clipping uses a boolean selection `(x >= lower) & (x < upper)`, which
  selects exactly what `searchsorted(side='left')` did on the sorted array.
  For float32 data, the limits are rounded up to the next float32
  (`_ceil_float32`), so the float32 comparisons give the same result as
  comparing in double precision, as the old `searchsorted` did.
- `d2` is recomputed from the current subset instead of being kept aligned
  with the sorted array.

Results (`compare_bstat.py`):

| | original | new |
| --- | --- | --- |
| per 150×150 box | 1.21 ms | 0.41 ms |
| full 5762² image with mask | 2.73 s | 1.15 s |

Iteration counts were identical for all sampled boxes. Differences come only
from the float32 summation order: median relative difference ~1e-7 to 1e-6;
max 6e-6 for the clipped rms. The clipped mean shows larger *relative*
differences (up to 1.7e-3) because its values are close to zero; relative to
the rms the difference is < 1e-6.

### 2.5 Results

`rms_mean_map` (`bench_rms_mean_map.py`, 20 cores, numpy 1.26.4):

| Case | Before | After |
| --- | --- | --- |
| 5762², no mask, box (150, 15) | 13.4 s | 5.7 s |
| 5762², masked region, box (60, 20) | 2.5 s | 1.4 s |
| 1200×1000 cutout, box (150, 50) (padded path) | 0.1 s | 0.1 s |
| 11524² (2×2 tiled), box (150, 15) | 60.7 s | 22.2 s |

Map differences (`compare_maps.py`): on the real image, all within 5e-5 of
the local rms, with identical non-finite entries. On the tiled image, 6 of
~580k boxes differ by more than 1e-4 relative. Five are ≤ 1e-3. One box next
to a bright source changed rms from 0.00051 to 0.00015: the clipping there
converges very slowly (65 iterations in the original; the new run reached
the fallback for > 198 iterations), so round-off decides the result. The
new value matches the neighbouring boxes (~0.00015).

End-to-end:

- `test/tbdsf_process_image.py`: passes with both HEAD and the new code.
  The two runs' outputs agree at rtol 1e-3, catalogs included.
- `benchmark/run_pybdsf.py` (NGC6048, adaptive rms box, wavelets):
  **2:20 → 1:38**. Rms map range (9.06e-05, 2.73e-03), 931 islands,
  26 small-scale sources and a total flux of 12.953 Jy are the same in both.

### 2.6 Not changed

- The per-island loops in `Op_rmsimage.__call__` (~10 s per call in the
  profile). Vectorizing them with `ndimage.maximum_position` would search
  the island pixels instead of the island's bounding box (as `nanargmax`
  over `image[s]` does now), which changes results.
- `multi_proc.parallel_map` itself still turns its input into an object
  array and uses a `Manager` queue. This is now cheap for `rms_mean_map`
  because the input holds only chunk numbers.

### 2.7 Files changed (uncommitted)

- `bdsf/functions.py`: `bstat` uses partition/boolean clipping; new helper
  `_ceil_float32`.
- `bdsf/rmsimage.py`: vectorized box setup, one `parallel_map` call per map,
  `process_mean_rms_boxes()`; removed `process_mean_rms_maps()`,
  `for_masked()` and `import itertools`.

---

## 3. Saving the work

**User:** Can you save the conversation above in a separate markdown file?
And can you also save benchmark scripts that you used for your benchmarks?
Lastly, can you also save additional tests you may have written?

**Claude:** Saved in `benchmark/rmsimage-optimization/` (see `README.md`).
No new unit tests were written; correctness was checked with the comparison
scripts and the existing end-to-end test, which `run_comparison.sh` reruns.

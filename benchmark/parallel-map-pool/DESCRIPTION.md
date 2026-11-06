# Rewrite of `multi_proc.parallel_map` with a process pool

Branch: `parallel-map-pool`, on top of `master` (5dab65e).

## Why `ProcessPoolExecutor` and not `Pool.imap`

The request was to use `multiprocessing.Pool.imap`. That has one problem: if
a `Pool` worker dies, for example because the out-of-memory killer stops it,
the pool starts a replacement and the dead worker's task is lost, so `imap`
waits forever. The old `parallel_map` checked worker exit codes and raised an
error. The rewrite therefore uses `concurrent.futures.ProcessPoolExecutor`.
It works the same way (a pool of forked workers that each take the next task
when free) but raises `BrokenProcessPool` when a worker dies. Switching to
`Pool.imap` is a small change if wanted.

## What changed

**Old behaviour:** `parallel_map` split the sequence into one fixed,
contiguous part per process, converted it to a numpy object array, and sent
results back through a `multiprocessing.Manager` queue.

**New behaviour:**
- **Dynamic scheduling:** tasks go to whichever worker is free.
  - With `weights` (gausfit), the most expensive items go first, one at a
    time.
  - Without weights, items go in batches, about 4 per worker.
- **Data:** the function and sequence are passed as per-pool initializer
  arguments, inherited through fork and stored in each worker's module-level
  variables. The caller's state is unchanged, so overlapping and nested calls
  remain independent. Nothing large is pickled, as before. Workers receive
  only item indices and send results back through the pool's pipes. The
  Manager queue and the object-array conversion are gone.
- **Status bar:** updated in the parent process as results arrive, instead of
  from the workers through a Manager dict.
- **Errors and interrupts:** if a worker raises an exception, a worker dies,
  or the call is interrupted (Ctrl+C), the pending tasks are cancelled and the
  remaining workers are terminated and joined before the exception is
  re-raised, as in the old implementation. `ProcessPoolExecutor.shutdown()`
  alone does not stop running tasks before Python 3.14, so this uses the
  executor's list of worker processes.
- **Callers:** gausfit, shapefit and psf_vary are unchanged. In rmsimage,
  only a comment in `process_mean_rms_boxes()` changed (see below).

`bdsf/multi_proc.py` is about 170 lines shorter.

## Testing

- **New test:** `test/tmulti_proc.py`. Its 12 checks all pass:
  - result order is kept, also for generator input and with weights
  - empty, single-item and two-item inputs
  - a closure over a large array works (inherited through fork, not pickled)
  - overlapping calls from threads and nested calls keep their task state separate
  - exceptions in a worker reach the caller
  - a killed worker raises `BrokenProcessPool` instead of hanging
  - an interrupted call (`KeyboardInterrupt`) leaves no workers running
  - the status bar is incremented once per item and stopped
  - the parent's module-level state remains unset afterwards

  The test was added to `.github/workflows/ci.yml` and to the cibuildwheel
  `test-command` in `pyproject.toml`.
- **End-to-end:** `test/tbdsf_process_image.py` passes. Its outputs match
  those of `master` at rtol 1e-3, catalogs included.
- **Benchmark:** `benchmark/run_pybdsf.py` (NGC6048, 20 cores, numpy 2.4.6,
  idle machine) took 1:04 on both `master` (5dab65e) and this branch. Both
  runs give the same rms map range (9.06e-05, 2.73e-03), 931 islands and a
  total flux of 12.908 Jy.
- **Chunk count in `process_mean_rms_boxes()`:** `master` splits the rms
  boxes into `100 * nproc` chunks, because the old `parallel_map` gave each
  process a fixed share. With dynamic scheduling fewer chunks might do, so
  `rms_mean_map` was timed on NGC6048 tiled 2x2
  (`benchmark/rmsimage-optimization/bench_rms_mean_map.py --tiled`,
  5 interleaved runs each):

  | Chunks per process | min | median |
  |---|---|---|
  | `master` (100) | 12.4 s | 14.0 s |
  | 100 | 12.8 s | 13.5 s |
  | 16 | 13.1 s | 14.4 s |
  | 4 | 13.7 s | 14.1 s |
  | 1 | 14.2 s | 15.1 s |

  Fewer chunks were not faster, so the count stays at 100; only the comment
  that explained it was updated. All maps were bit-identical to `master`.

On this machine, the multiprocessing overhead is small compared with the
work done in the workers, so the rewrite gives no measurable speedup. A gain
on much larger images (such as in the rapthor profile,
`benchmark/filter_skymodel.prof`) has not been measured. The main benefits
are simpler code, dynamic load balancing, and an error instead of a hang
when a worker dies.

## Notes

- Rebased onto `master` after #504 (`rmsimage-optimization`) was merged.
- To reproduce the checks: `benchmark/parallel-map-pool/run_comparison.sh`
  (baseline: `master`). `setup_overlays.sh` now also copies `bdsf.libs`
  (the libgfortran bundled with the wheel).

# Rewrite of `multi_proc.parallel_map` with a process pool

Branch: `parallel-map-pool`, on top of `master`.

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
- **Data:** the function and sequence are stored in module-level variables
  before the workers are forked, so the workers inherit them and nothing
  large is pickled, as before. Workers receive only item indices and send
  results back through the pool's pipes. The Manager queue and the
  object-array conversion are gone.
- **Status bar:** updated in the parent process as results arrive, instead of
  from the workers through a Manager dict.
- **Errors and interrupts:** if a worker raises an exception, a worker dies,
  or the call is interrupted (Ctrl+C), the pending tasks are cancelled and the
  remaining workers are terminated and joined before the exception is
  re-raised, as in the old implementation. `ProcessPoolExecutor.shutdown()`
  alone does not stop running tasks before Python 3.14, so this uses the
  executor's list of worker processes.
- **Callers:** gausfit, shapefit, psf_vary and rmsimage are unchanged.

`bdsf/multi_proc.py` is about 170 lines shorter.

## Testing

- **New test:** `test/tmulti_proc.py`. Its 10 checks all pass:
  - result order is kept, also for generator input and with weights
  - empty, single-item and two-item inputs
  - a closure over a large array works (inherited through fork, not pickled)
  - exceptions in a worker reach the caller
  - a killed worker raises `BrokenProcessPool` instead of hanging
  - an interrupted call (`KeyboardInterrupt`) leaves no workers running
  - the status bar is incremented once per item and stopped
  - the module-level state is cleared afterwards

  The test was added to `.github/workflows/ci.yml` and to the cibuildwheel
  `test-command` in `pyproject.toml`.
- **End-to-end:** `test/tbdsf_process_image.py` passes. Its outputs match
  those of `master` at rtol 1e-3, catalogs included.
- **Benchmark:** `benchmark/run_pybdsf.py` (NGC6048, 20 cores) took 2:17 on
  `master` and 2:13 with this branch, within noise (the machine had a load
  average of 5-7 from other users during the runs). Both runs give the same
  rms map range (9.06e-05, 2.73e-03), 931 islands and a total flux of
  12.953 Jy.

On this machine, the multiprocessing overhead is small compared with the
work done in the workers. The gain should show on the much larger images
seen in the rapthor profile (`benchmark/filter_skymodel.prof`), but that has
not been measured. The main benefits are simpler code, dynamic load
balancing, and an error instead of a hang when a worker dies.

## Notes

- This branch is independent of `rmsimage-optimization`; the two can be
  merged in either order. `rms_mean_map` on `master` still sends one item per
  box, which the new `parallel_map` hands out in batches.
- To reproduce the checks: `benchmark/parallel-map-pool/run_comparison.sh`
  (baseline: `master`).

"""Multiprocessing module to handle parallelization.

This module can optionally update a statusbar and can use weights to
schedule the most expensive tasks first.

The tasks are distributed dynamically over a pool of worker processes, so
that a worker that finishes early picks up the next task. The worker
processes inherit the function and sequence through fork, and a per-pool
initializer stores them in each worker: only item indices are sent to the
workers and only the results are sent back.
"""

import concurrent.futures
import concurrent.futures.process
import math
import multiprocessing
import os

import numpy

def nproc():
    """
    Return the number of CPU cores _available_ to the current process, similar
    to what the Linux `nproc` command does. This can be less than the total
    number of CPU cores in the machine.
    NOTE: This function uses `os.sched_getaffinity()`, which is not available
    on every OS. Use `os.cpu_count()` as fall-back; return 1 in the rare case
    that `os.cpu_count()` returns `None`.
    """
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1

_ncpus = nproc()

# PyBDSF currently relies on fork-style multiprocessing.
# The spawn and forkserver start methods are not supported because
# parts of the codebase assume inherited interpreter state and are
# not safe for re-import during process startup.
fork_context = multiprocessing.get_context("fork")

__all__ = ('parallel_map',)

# Task state is assigned only in workers, never in the calling process.
# This keeps concurrent pools and calls nested inside a worker independent.
_task_function = None
_task_sequence = None


def _init_worker(function, sequence):
    """Set this pool's task state using arguments inherited through fork."""
    global _task_function, _task_sequence
    _task_function, _task_sequence = function, sequence
    _forget_inherited_executors()


def _forget_inherited_executors():
    """
    Forget the process pools that this (forked) worker inherited from its
    parent.

    `concurrent.futures.process` keeps a registry of the "wakeup" objects of
    all running process pools, and a forked worker inherits it. When the
    worker exits, it wakes all of them up, which takes the lock of each
    wakeup object. If another thread of the parent held one of those locks
    when the worker was forked (e.g., when parallel_map is called from two
    threads at once), the worker inherits a lock that is never released, and
    hangs when it exits; so does parallel_map, while waiting for it. This
    affects recent Python versions, in which these objects have a lock (e.g.,
    3.12.15, 3.13 and 3.14, but not 3.12.3).
    The registry of the parent's pools is of no use in a worker, so clear it,
    like CPython does for thread pools in forked processes.
    """
    wakeups = getattr(concurrent.futures.process, '_threads_wakeups', None)
    if wakeups is not None:
        wakeups.clear()


def _run_batch(indices):
    """
    Applies the current task function to the items of the current task
    sequence with the given indices (in a worker process).

    :param indices: indices of the items to process
    :returns: list of (index, result) tuples
    """
    return [(idx, _task_function(_task_sequence[idx])) for idx in indices]


def _terminate_workers(executor):
    """Cancel the pending tasks of the executor and stop its workers"""
    # ProcessPoolExecutor has no public way to stop running tasks before
    # Python 3.14 (terminate_workers()), so use its list of worker processes
    processes = list((executor._processes or {}).values())
    executor.shutdown(wait=False, cancel_futures=True)
    for process in processes:
        if process.is_alive():
            process.terminate()
    for process in processes:
        process.join()


def parallel_map(function, sequence, numcores=None, bar=None, weights=None):
    """
    A parallelized version of the native Python map function that
    utilizes the Python multiprocessing module to divide and
    conquer a sequence.

    parallel_map does not yet support multiple argument sequences.

    :param function: callable function that accepts argument from iterable
    :param sequence: iterable sequence
    :param numcores: number of cores to use (if None, all are used)
    :param bar: statusbar to update during fit
    :param weights: estimated cost of each item; items are processed in order
        of decreasing weight, one at a time, which keeps the workers evenly
        loaded when the costs vary a lot

    """
    if not callable(function):
        raise TypeError("input function '%s' is not callable" %
                        repr(function))

    if not numpy.iterable(sequence):
        raise TypeError("input '%s' is not iterable" %
                        repr(sequence))

    sequence = list(sequence)
    size = len(sequence)

    if numcores is None:
        numcores = _ncpus - 1
    numcores = max(1, min(numcores, _ncpus - 1, size))

    # Divide the items into batches, which are handed out to the workers as
    # they become free. With weights, process items one at a time, most
    # expensive first; otherwise, use a few batches per worker to limit the
    # communication overhead.
    if weights is None:
        order = numpy.arange(size)
        batch_size = max(1, math.ceil(size / (4 * numcores)))
    else:
        order = numpy.argsort(-numpy.asarray(weights, dtype=numpy.float64),
                              kind='stable')
        batch_size = 1
    batches = [order[start:start + batch_size].tolist()
               for start in range(0, size, batch_size)]

    results = [None] * size
    executor = concurrent.futures.ProcessPoolExecutor(max_workers=numcores,
                                                      mp_context=fork_context,
                                                      initializer=_init_worker,
                                                      initargs=(function, sequence))
    try:
        futures = [executor.submit(_run_batch, batch) for batch in batches]
        for future in concurrent.futures.as_completed(futures):
            for idx, result in future.result():
                results[idx] = result
                if bar is not None and bar.started:
                    bar.increment()
    except BaseException:
        # An exception (e.g., KeyboardInterrupt, or BrokenProcessPool after a
        # worker was killed by the out-of-memory killer) would otherwise leave
        # the remaining workers running until their tasks finish
        _terminate_workers(executor)
        raise
    executor.shutdown()

    if bar is not None and bar.started:
        bar.stop()
    return results

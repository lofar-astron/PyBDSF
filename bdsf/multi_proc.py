"""Multiprocessing module to handle parallelization.

This module can optionally update a statusbar and can divide tasks
between cores using weights (so that each core gets a set of tasks with
the same total weight).

Adapted from a module by Brian Refsdal at SAO, available at AstroPython
(http://www.astropython.org/snippet/2010/3/Parallel-map-using-multiprocessing).

"""

import heapq
import multiprocessing
import os
import sys
import traceback

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


def worker(f, ii, chunk, out_q, err_q, lock, bar, bar_state, preserve_order=False):
    """
    A worker function that maps an input function over a
    slice of the input iterable.

    :param f  : callable function that accepts argument from iterable
    :param ii  : process ID
    :param chunk: slice of input iterable
    :param out_q: process-safe output queue
    :param err_q: process-safe queue to populate on exception
    :param lock : lock object (kept for signature compatibility)
    :param bar: statusbar to update during fit
    :param bar_state: dictionary holding shared memory Values for statusbar state
    :param preserve_order: whether chunk entries carry their original index
    """
    vals = []

    for entry in chunk:
        if preserve_order:
            val_idx, val = entry
        else:
            val = entry
        try:
            result = f(val)
        except Exception as e:
            print("Thread raised exception", e, file=sys.stderr)
            print("Traceback of thread is:", file=sys.stderr)
            print("-------------------------", file=sys.stderr)
            traceback.print_exception(e, file=sys.stderr)
            print("-------------------------", file=sys.stderr)
            err_q.put(e)
            return

        if preserve_order:
            vals.append((val_idx, result))
        else:
            vals.append(result)

        if bar is not None:
            # Access shared memory values directly without manager overhead
            if bool(bar_state['started'].value):
                bar.pos = bar_state['pos'].value
                bar.spin_pos = bar_state['spin_pos'].value
                bar.started = bool(bar_state['started'].value)
                increment = bar.increment()
                bar_state['started'].value = int(bar.started)
                bar_state['pos'].value += increment
                bar_state['spin_pos'].value += increment
                if bar_state['spin_pos'].value >= 4:
                    bar_state['spin_pos'].value = 0

    out_q.put((ii, vals))


def run_tasks(procs, err_q, out_q, num, preserve_order=False, total_items=None):
    """
    A function that executes populated processes and processes
    the resultant array. Checks error queue for any exceptions.

    :param procs: list of Process objects
    :param err_q: thread-safe queue to populate on exception
    :param out_q: thread-safe output queue
    :param num : length of resultant array
    :param preserve_order: whether worker outputs carry original item indices
    :param total_items: total number of items to reconstruct

    """
    die = (lambda vals: [val.terminate() for val in vals
                         if val.exitcode is None])

    try:
        for proc in procs:
            proc.start()

        # Retrieve results from queue BEFORE joining processes to avoid
        # IPC pipe buffer deadlocks with context.Queue
        raw_results = []
        for _ in range(num):
            raw_results.append(out_q.get())

        for proc in procs:
            proc.join()
            if proc.exitcode != 0:
                raise RuntimeError(
                    f"Process {proc.name} died unexpectedly (exit code: {proc.exitcode})"
                )

    except Exception as e:
        die(procs)
        raise e

    if not err_q.empty():
        die(procs)
        raise err_q.get()

    if preserve_order:
        results = [None] * total_items
        for idx, result in raw_results:
            for item_idx, item_result in result:
                results[item_idx] = item_result
        return results

    results = [None] * num
    for idx, result in raw_results:
        results[idx] = result

    result_list = []
    for result in results:
        result_list += result

    return result_list


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
    :param weights: weights to use when splitting the sequence

    """
    if not callable(function):
        raise TypeError("input function '%s' is not callable" %
                        repr(function))

    if not numpy.iterable(sequence):
        raise TypeError("input '%s' is not iterable" %
                        repr(sequence))

    # Removed numpy.array(..., dtype=object) casting to prevent deep inspection overhead
    sequence_list = list(sequence)
    size = len(sequence_list)

    if size == 1:
        results = list(map(function, sequence_list))
        if bar is not None:
            bar.stop()
        return results

    if numcores is None:
        numcores = _ncpus - 1
    if numcores > _ncpus - 1:
        numcores = _ncpus - 1
    if numcores < 1:
        numcores = 1

    # Use fast context-native queues instead of Manager proxy queues
    out_q = fork_context.Queue()
    err_q = fork_context.Queue()
    lock = None

    # Use un-locked shared memory Value objects for statusbar state
    # to eliminate IPC Manager overhead
    bar_state = {}
    if bar is not None:
        bar_state['pos'] = fork_context.Value('i', int(bar.pos), lock=False)
        bar_state['spin_pos'] = fork_context.Value('i', int(bar.spin_pos), lock=False)
        bar_state['started'] = fork_context.Value('i', int(bar.started), lock=False)

    if size < numcores:
        numcores = size

    preserve_order = False
    if weights is None or numcores == size:
        # Native Python split equivalent to numpy.array_split
        # (no numpy overhead)
        n_each_section, extras = divmod(size, numcores)
        section_sizes = ([n_each_section + 1] * extras) + ([n_each_section] * (numcores - extras))
        
        sequence = []
        start = 0
        for sec_size in section_sizes:
            sequence.append(sequence_list[start:start + sec_size])
            start += sec_size
    else:
        preserve_order = True
        weight_array = numpy.asarray(weights, dtype=numpy.float64)
        indexed_sequence = list(enumerate(sequence_list))
        bins = [[] for _ in range(numcores)]
        heap = [(0.0, idx) for idx in range(numcores)]
        heapq.heapify(heap)

        weighted_items = zip(indexed_sequence, weight_array.tolist())
        for (orig_idx, item), weight in sorted(weighted_items,
                                              key=lambda pair: pair[1],
                                              reverse=True):
            current_weight, worker_idx = heapq.heappop(heap)
            bins[worker_idx].append((orig_idx, item))
            current_weight += float(weight)
            heapq.heappush(heap, (current_weight, worker_idx))

        sequence = bins

    while len(sequence[-1]) == 0:
        sequence.pop()

    procs = [fork_context.Process(target=worker,
             args=(function, ii, chunk, out_q, err_q, lock, bar, bar_state,
                   preserve_order))
             for ii, chunk in enumerate(sequence)]

    try:
        results = run_tasks(procs, err_q, out_q, len(sequence),
                            preserve_order=preserve_order, total_items=size)
        if bar is not None:
            if bar.started:
                bar.stop()
        return results

    except KeyboardInterrupt:
        for proc in procs:
            if proc.exitcode is None:
                proc.terminate()
                proc.join()
        raise
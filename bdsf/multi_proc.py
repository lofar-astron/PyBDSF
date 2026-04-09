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

# Try to determine the number of CPU cores _available_ to the current process,
# similar to what the Linux `nproc` command does. If that fails, return the
# total number of CPU cores in the machine.
try:
    _ncpus = len(os.sched_getaffinity(0))
except AttributeError:
    _ncpus = multiprocessing.cpu_count()

# Set the start method to "fork". Other methods don't work with our codebase.
multiprocessing.set_start_method('fork')

__all__ = ('parallel_map',)


def worker(f, ii, chunk, out_q, err_q, lock, bar, bar_state, preserve_order=False):
    """
    A worker function that maps an input function over a
    slice of the input iterable.

    :param f  : callable function that accepts argument from iterable
    :param ii  : process ID
    :param chunk: slice of input iterable
    :param out_q: thread-safe output queue
    :param err_q: thread-safe queue to populate on exception
    :param lock : thread-safe lock to protect a resource
           ( useful in extending parallel_map() )
    :param bar: statusbar to update during fit
    :param bar_state: statusbar state dictionary
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
            etype, val, tbk = sys.exc_info()
            print('Thread raised exception', e)
            print('Traceback of thread is:')
            print('-------------------------')
            traceback.print_tb(tbk)
            print('-------------------------')
            err_q.put(e)
            return

        if preserve_order:
            vals.append((val_idx, result))
        else:
            vals.append(result)

        if bar is not None:
            if bar_state['started']:
                bar.pos = bar_state['pos']
                bar.spin_pos = bar_state['spin_pos']
                bar.started = bar_state['started']
                increment = bar.increment()
                bar_state['started'] = bar.started
                bar_state['pos'] += increment
                bar_state['spin_pos'] += increment
                if bar_state['spin_pos'] >= 4:
                    bar_state['spin_pos'] = 0

    out_q.put((ii, vals))


def run_tasks(procs, err_q, out_q, num, preserve_order=False, total_items=None):
    """
    A function that executes populated processes and processes
    the resultant array. Checks error queue for any exceptions.

    :param procs: list of Process objects
    :param out_q: thread-safe output queue
    :param err_q: thread-safe queue to populate on exception
    :param num : length of resultant array
    :param preserve_order: whether worker outputs carry original item indices
    :param total_items: total number of items to reconstruct

    """
    die = (lambda vals: [val.terminate() for val in vals
                         if val.exitcode is None])

    try:
        for proc in procs:
            proc.start()

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
        for i in range(num):
            idx, result = out_q.get()
            for item_idx, item_result in result:
                results[item_idx] = item_result
        return results

    results = [None] * num
    for i in range(num):
        idx, result = out_q.get()
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

    sequence_list = list(sequence)
    sequence = numpy.array(sequence_list, dtype=object)
    size = len(sequence)

    if size == 1:
        results = list(map(function, sequence))
        if bar is not None:
            bar.stop()
        return results

    if numcores is None:
        numcores = _ncpus - 1
    if numcores > _ncpus - 1:
        numcores = _ncpus - 1
    if numcores < 1:
        numcores = 1

    manager = multiprocessing.Manager()

    out_q = manager.Queue()
    err_q = manager.Queue()
    lock = manager.Lock()
    bar_state = manager.dict()
    if bar is not None:
        bar_state['pos'] = bar.pos
        bar_state['spin_pos'] = bar.spin_pos
        bar_state['started'] = bar.started

    if size < numcores:
        numcores = size

    preserve_order = False
    if weights is None or numcores == size:
        sequence = numpy.array_split(sequence, numcores)
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

    procs = [multiprocessing.Process(target=worker,
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

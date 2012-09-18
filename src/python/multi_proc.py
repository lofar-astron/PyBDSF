"""Multiprocessing module to handle parallelization.

This module can optionally update a statusbar and can divide tasks
between cores using weights (so that each core gets a set of tasks with
the same total weight).

Adapted from a module by Brian Refsdal at SAO, available at AstroPython
(http://www.astropython.org/snippet/2010/3/Parallel-map-using-multiprocessing).

"""
import numpy
_multi=False
_ncpus=1

try:
    # May raise ImportError
    import multiprocessing
    _multi=True

    # May raise NotImplementedError
    _ncpus = multiprocessing.cpu_count()
except:
    pass


__all__ = ('parallel_map',)


def worker(f, ii, chunk, out_q, err_q, lock, bar, bar_state):
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
    """
    vals = []

    # iterate over slice
    for val in chunk:
        try:
            result = f(val)
        except Exception, e:
            err_q.put(e)
            return

        vals.append(result)

        # update statusbar
        if bar != None:
            if bar_state['started']:
                bar.pos = bar_state['pos']
                bar.spin_pos = bar_state['spin_pos']
                bar.increment()
                bar_state['pos'] += 1
                bar_state['spin_pos'] += 1
                if bar_state['spin_pos'] >= 4:
                    bar_state['spin_pos'] = 0

    # output the result and task ID to output queue
    out_q.put( (ii, vals) )


def run_tasks(procs, err_q, out_q, num):
    """
    A function that executes populated processes and processes
    the resultant array. Checks error queue for any exceptions.

    :param procs: list of Process objects
    :param out_q: thread-safe output queue
    :param err_q: thread-safe queue to populate on exception
    :param num : length of resultant array

    """
    # function to terminate processes that are still running.
    die = (lambda vals : [val.terminate() for val in vals
               if val.exitcode is None])

    try:
        for proc in procs:
            proc.start()

        for proc in procs:
            proc.join()

    except Exception, e:
        # kill all slave processes on ctrl-C
        die(procs)
        raise e

    if not err_q.empty():
        # kill all on any exception from any one slave
        die(procs)
        raise err_q.get()

    # Processes finish in arbitrary order. Process IDs double
    # as index in the resultant array.
    results=[None]*num;
    for i in range(num):
        idx, result = out_q.get()
        results[idx] = numpy.array(result, dtype=object)

    # Remove extra dimension added by array_split
    return numpy.concatenate(results).tolist()



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

    sequence = list(sequence)
    size = len(sequence)

    if not _multi or size == 1:
        return map(function, sequence)

    # Set default number of cores to use. Leave one core free for pyplot.
    if numcores is None:
        numcores = _ncpus - 1
    if numcores > _ncpus - 1:
        numcores = _ncpus - 1

    # Returns a started SyncManager object which can be used for sharing
    # objects between processes. The returned manager object corresponds
    # to a spawned child process and has methods which will create shared
    # objects and return corresponding proxies.
    manager = multiprocessing.Manager()

    # Create FIFO queue and lock shared objects and return proxies to them.
    # The managers handles a server process that manages shared objects that
    # each slave process has access to. Bottom line -- thread-safe.
    out_q = manager.Queue()
    err_q = manager.Queue()
    lock = manager.Lock()
    bar_state = manager.dict()
    if bar != None:
        bar_state['pos'] = bar.pos
        bar_state['spin_pos'] = bar.spin_pos
        bar_state['started'] = bar.started

    # if sequence is less than numcores, only use len sequence number of
    # processes
    if size < numcores:
        numcores = size

    # group sequence into numcores-worth of chunks
    if weights == None or numcores == size:
        # No grouping specified (or there are as many cores as
        # processes), so divide into equal chunks
        sequence = numpy.array_split(sequence, numcores)
    else:
        # Group so that each group has roughly an equal sum of weights
        weight_per_core = numpy.sum(weights)/float(numcores)
        cut_values = []
        temp_sum = 0.0
        for indx, weight in enumerate(weights):
            temp_sum += weight
            if temp_sum > weight_per_core:
                cut_values.append(indx)
                temp_sum = weight
        if len(cut_values) > numcores - 1:
            cut_values = cut_values[0:numcores-1]
        sequence = numpy.array_split(sequence, cut_values)

    procs = [multiprocessing.Process(target=worker,
             args=(function, ii, chunk, out_q, err_q, lock, bar, bar_state))
             for ii, chunk in enumerate(sequence)]

    try:
        results = run_tasks(procs, err_q, out_q, len(sequence))
        if bar != None:
            if bar.started:
                bar.pos = bar_state['pos']
                bar.spin_pos = bar_state['spin_pos']
                while bar.pos < bar.max:
                    bar.increment()
        return results

    except KeyboardInterrupt:
        for proc in procs:
            if proc.exitcode is None:
                proc.terminate()
                proc.join()
        raise

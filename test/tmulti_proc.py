"""Tests for bdsf.multi_proc.parallel_map

Run with: python tmulti_proc.py (exit status 0 if all tests pass)
"""
import concurrent.futures
import multiprocessing
import os
import signal
import subprocess
import sys
import threading
import time
from unittest.mock import patch

import numpy as np

from bdsf import multi_proc as mp


class FakeBar:
    """Minimal stand-in for bdsf.statusbar.StatusBar"""
    def __init__(self):
        self.started = 1
        self.pos = 0
        self.stopped = False

    def increment(self):
        self.pos += 1

    def stop(self):
        self.stopped = True
        self.started = 0


def square(x):
    return x * x


def fail_on_five(x):
    if x == 5:
        raise ValueError('item 5 failed')
    return x


def die_on_three(x):
    if x == 3:
        os._exit(1)  # simulate a worker killed by, e.g., the OOM killer
    time.sleep(0.05)
    return x


def test_results_in_order():
    items = list(range(1000))
    assert mp.parallel_map(square, items, numcores=4) == [x * x for x in items]


def test_generator_input():
    assert mp.parallel_map(square, (x for x in range(10)), numcores=3) == \
        [x * x for x in range(10)]


def test_weights_keep_order():
    items = list(range(50))
    weights = np.random.default_rng(1).random(50)
    assert mp.parallel_map(square, items, numcores=4, weights=weights) == \
        [x * x for x in items]


def test_small_sequences():
    assert mp.parallel_map(square, [], numcores=4) == []
    assert mp.parallel_map(square, [3], numcores=4) == [9]
    assert mp.parallel_map(square, [1, 2], numcores=8) == [1, 4]


def test_closure_inherits_data():
    # The workers are forked, so a closure over large data works without
    # pickling it
    data = np.arange(10_000_000, dtype=np.float64)
    result = mp.parallel_map(lambda k: data[k::4].sum(), range(4), numcores=4)
    assert result == [data[k::4].sum() for k in range(4)]


def test_overlapping_calls():
    # Hold both callers before their first submission so their pools start
    # after both calls have prepared their functions and input sequences.
    ready = threading.Barrier(2, timeout=10)

    class SynchronizedExecutor(concurrent.futures.ProcessPoolExecutor):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            ready.wait()

    def run(offset):
        return mp.parallel_map(lambda x: x + offset,
                               range(offset, offset + 8), numcores=2)

    with patch.object(concurrent.futures, 'ProcessPoolExecutor',
                      SynchronizedExecutor):
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            calls = [executor.submit(run, offset) for offset in (100, 200)]
            results = [call.result() for call in calls]
    assert results == [list(range(200, 208)), list(range(400, 408))]


# Run in a separate process, because it hangs if the bug is present
FORK_WHILE_LOCKED = """
import concurrent.futures
from bdsf import multi_proc as mp
other = concurrent.futures.ProcessPoolExecutor(1, mp_context=mp.fork_context)
other.submit(int).result()  # starts its manager thread
wakeup = other._executor_manager_thread_wakeup
if not hasattr(wakeup, '_lock'):  # Python versions without the lock
    print('skipped')
else:
    # As if another thread is waking up `other` while parallel_map forks
    with wakeup._lock:
        result = mp.parallel_map(abs, range(-4, 4), numcores=2)
    assert result == [4, 3, 2, 1, 0, 1, 2, 3], result
    print('ok')
other.shutdown()
"""


def test_fork_while_other_pool_is_busy():
    # Forked workers inherit the wakeup objects of all process pools of the
    # parent, and wake them up when they exit. If another thread held the
    # lock of such an object during the fork (which test_overlapping_calls
    # can trigger), the workers hung when exiting, and so did parallel_map.
    proc = subprocess.Popen([sys.executable, '-c', FORK_WHILE_LOCKED],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, start_new_session=True)
    try:
        output, _ = proc.communicate(timeout=60)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)  # including the hung workers
        proc.communicate()
        raise AssertionError('parallel_map hung')
    assert proc.returncode == 0, output
    assert output.splitlines()[-1] in ('ok', 'skipped'), output


def test_nested_calls():
    def outer(x):
        return sum(mp.parallel_map(lambda y: y + x, [1, 2], numcores=1))

    # A single outer worker must retain its state across items and batches.
    assert mp.parallel_map(outer, range(8), numcores=1) == \
        [3 + 2 * x for x in range(8)]


def test_exception_propagates():
    try:
        mp.parallel_map(fail_on_five, range(20), numcores=4)
    except ValueError as e:
        assert 'item 5 failed' in str(e)
    else:
        raise AssertionError('ValueError not raised')


def test_killed_worker_raises():
    start = time.time()
    try:
        mp.parallel_map(die_on_three, range(40), numcores=4)
    except concurrent.futures.process.BrokenProcessPool:
        pass
    else:
        raise AssertionError('BrokenProcessPool not raised')
    assert time.time() - start < 30


def test_interrupt_stops_workers():
    def interrupt(signum, frame):
        raise KeyboardInterrupt
    old_handler = signal.signal(signal.SIGALRM, interrupt)
    signal.alarm(2)
    start = time.time()
    try:
        mp.parallel_map(time.sleep, [60] * 4, numcores=2)
    except KeyboardInterrupt:
        pass
    else:
        raise AssertionError('KeyboardInterrupt not raised')
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)
    assert time.time() - start < 30
    assert not multiprocessing.active_children(), 'workers still running'


def test_statusbar():
    bar = FakeBar()
    mp.parallel_map(square, range(100), numcores=4, bar=bar)
    assert bar.pos == 100
    assert bar.stopped


if __name__ == '__main__':
    failed = 0
    for name, test in list(globals().items()):
        if name.startswith('test_') and callable(test):
            try:
                test()
                print('PASSED', name)
            except Exception as e:
                failed += 1
                print('FAILED', name, repr(e))
    sys.exit(1 if failed else 0)

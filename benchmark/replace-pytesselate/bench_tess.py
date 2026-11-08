"""Compare the results, speed and memory use of `_pytesselate` and `tessellate.py`.

Usage: bench_tess.py <tessellate.py> [records.pkl ...] [--repeat N] [--synthetic]

For each recorded call (see `record_tess.py`), the generators are tessellated
with all methods that PyBDSF can use (`pytess_simple` with the weights
'unity', 'log10' and 'sqrtlog10', and `pytess_roundness`), each with code 's'
(normal) and 'c' (fuzzy). With --synthetic, cases with random generators
(more than PyBDSF typically uses) are added.

Each case is run in a separate process for each implementation, which
reports the best time of at least N calls (default 3; more for fast cases, so
that they take at least 1 s in total), and the increase of the peak resident
set size (VmHWM) during the first call. The results of both implementations must
be identical. Fuzzy `pytess_roundness` is not compared, because its result is
undefined in the Fortran code (see `tessellate.py`).
"""
import importlib.util
import json
import os
import pickle
import subprocess
import sys
import tempfile
import time

import numpy as N


def load_python(path):
    spec = importlib.util.spec_from_file_location('tessellate', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def memory(field):
    """VmRSS (current) or VmHWM (peak) resident set size of this process, in bytes."""
    with open('/proc/self/status') as f:
        for line in f:
            if line.startswith(field + ':'):
                return int(line.split()[1]) * 1024


def worker(impl, pyfile, case_json, out, repeat):
    """Run one case with one implementation (in a separate process)."""
    case = json.loads(case_json)
    # bdsf is imported in both cases, as in a PyBDSF run. This matters for
    # small images: in a process with a small heap, glibc returns the memory
    # freed after each call to the system, and the next call takes page
    # faults (which makes the Python version up to 30% slower).
    import bdsf
    if impl == 'fortran':
        from bdsf import _pytesselate as module
    else:
        module = load_python(pyfile)
    func = getattr(module, case['func'])
    args = case['args']
    # Touch everything that the call needs, except its own allocations
    N.zeros(10)
    # Reset the peak RSS (VmHWM) to the current RSS
    with open('/proc/self/clear_refs', 'w') as f:
        f.write('5')
    rss0 = memory('VmRSS')
    times = []
    t0 = time.perf_counter()
    result = func(*args)
    times.append(time.perf_counter() - t0)
    peak = memory('VmHWM') - rss0
    N.save(out, result)
    del result
    while len(times) < repeat or sum(times) < 1.0:
        t0 = time.perf_counter()
        func(*args)
        times.append(time.perf_counter() - t0)
    print(json.dumps(dict(time=min(times), peak=peak)))


def _fmt(t):
    return f'{t:7.3f}s' if t >= 0.1 else f'{t * 1e3:6.2f}ms'


def make_cases(records, synthetic):
    cases = []
    for name, rec in records:
        n, m, xg, yg, snr = rec['args'][:5]
        xg, yg, snr = (list(map(float, a)) for a in (xg, yg, snr))
        weights = {'unity': N.ones(len(snr)), 'log10': N.log10(snr),
                   'sqrtlog10': N.sqrt(N.log10(snr))}
        for code in 'sc':
            for wname, w in weights.items():
                cases.append((f'{name} simple {wname} {code}', 'pytess_simple',
                              [n, m, xg, yg, snr, list(map(float, w)), 0.05, code]))
            cases.append((f'{name} roundness {code}', 'pytess_roundness',
                          [n, m, xg, yg, snr, 0.05, code]))
    if synthetic:
        rng = N.random.default_rng(1)
        for n, ngens in ((5762, 50), (5762, 200), (2000, 1000)):
            xg = list(rng.uniform(0, n, ngens))
            yg = list(rng.uniform(0, n, ngens))
            snr = list(rng.uniform(10, 1000, ngens))
            w = list(N.log10(snr))
            name = f'random {n}x{n}, {ngens} gens'
            for code in 'sc':
                cases.append((f'{name} simple log10 {code}', 'pytess_simple',
                              [n, n, xg, yg, snr, w, 0.05, code]))
            cases.append((f'{name} roundness s', 'pytess_roundness',
                          [n, n, xg, yg, snr, 0.05, 's']))
            # A weight of 0 (invalid), for which the distances cannot be
            # bounded
            w0 = [0.0] + w[1:]
            cases.append((f'{name} simple w[0]=0 s', 'pytess_simple',
                          [n, n, xg, yg, snr, w0, 0.05, 's']))
    return cases


def main():
    args = sys.argv[1:]
    repeat = 3
    if '--repeat' in args:
        i = args.index('--repeat')
        repeat = int(args[i + 1])
        del args[i:i + 2]
    synthetic = '--synthetic' in args
    args = [a for a in args if a != '--synthetic']
    pyfile, recfiles = args[0], args[1:]

    records = []
    for recfile in recfiles:
        with open(recfile, 'rb') as f:
            for i, rec in enumerate(pickle.load(f)):
                records.append((f'{os.path.basename(recfile)}[{i}]', rec))
    cases = make_cases(records, synthetic)

    print(f'{"case":<42} {"shape":>11} {"ngens":>5}  {"Fortran":>8} {"Python":>8} '
          f'{"ratio":>6}  {"peak MB F":>9} {"peak MB P":>9}  identical')
    all_ok = True
    with tempfile.TemporaryDirectory() as tmp:
        for name, func, cargs in cases:
            res = {}
            for impl in ('fortran', 'python'):
                out = os.path.join(tmp, f'{impl}.npy')
                proc = subprocess.run(
                    [sys.executable, '-I', __file__, '--worker', impl, pyfile,
                     json.dumps(dict(func=func, args=cargs)), out, str(repeat)],
                    capture_output=True, text=True)
                if proc.returncode != 0:
                    res[impl] = None
                    print(f'{name}: {impl} failed: {proc.stderr.strip()[-300:]}')
                    continue
                res[impl] = json.loads(proc.stdout.strip().splitlines()[-1])
                res[impl]['out'] = N.load(out)
            undefined = func == 'pytess_roundness' and cargs[-1] != 's'
            if res['fortran'] is None and undefined:
                print(f'{name}: (expected, since the Fortran result is undefined)')
                continue
            if res['fortran'] is None or res['python'] is None:
                all_ok = False
                continue
            if undefined:
                same = 'n/a (undefined)'
            else:
                ok = N.array_equal(res['fortran']['out'], res['python']['out'])
                all_ok &= ok
                ndiff = N.count_nonzero(res['fortran']['out'] != res['python']['out'])
                same = 'yes' if ok else f'NO ({ndiff} pixels differ)'
            tf, tp = res['fortran']['time'], res['python']['time']
            print(f'{name:<42} {cargs[0]:>5}x{cargs[1]:<5} {len(cargs[2]):>5}  '
                  f'{_fmt(tf)} {_fmt(tp)} {tf / tp:6.1f}  '
                  f'{res["fortran"]["peak"] / 2**20:9.0f} {res["python"]["peak"] / 2**20:9.0f}  {same}',
                  flush=True)
    print('All results identical' if all_ok else 'DIFFERENCES FOUND')
    return 0 if all_ok else 1


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--worker':
        worker(*sys.argv[2:6], int(sys.argv[6]))
    else:
        sys.exit(main())

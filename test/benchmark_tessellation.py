"""Compare tessellation runtimes after checking exact output equivalence.

Run: python test/benchmark_tessellation.py --output tessellation-benchmark.json
Compilation, imports, data generation, equality checks and warmup are not timed.
Both implementations use the same inputs; output allocation is included.
"""

import argparse
import hashlib
import json
import platform
from statistics import median
import subprocess
import sys
from time import perf_counter

import numpy as np

from tessellation_reference import ROOT, implementations


def main():
    """Benchmark the requested size/count matrix and optionally save JSON.

    Speedup is median Fortran time divided by median Python time, so values
    above one favor Python. Every timed workload must first pass exact parity.
    ``--min-speedup`` optionally enforces a per-workload performance floor;
    raw samples are written before reporting a missed floor.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=[512, 1024])
    parser.add_argument("--generators", type=int, nargs="+", default=[2, 8, 64, 256])
    parser.add_argument("--repeat", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260921)
    parser.add_argument("--output")
    parser.add_argument("--min-speedup", type=float, help="Fail if any median speedup is below this value")
    options = parser.parse_args()
    if options.repeat < 1 or min(options.sizes + options.generators) < 1:
        parser.error("Sizes, generator counts and repeat must be positive")
    python, fortran = implementations()
    rng = np.random.default_rng(options.seed)
    rows = []
    print("shape       generators mode        Fortran ms   Python ms   speedup", flush=True)
    for size in options.sizes:
        for count in options.generators:
            x, y = rng.uniform(1, size, (2, count))
            snr = np.full(count, 10.)
            weights = rng.uniform(0.5, 2., count)
            # Reuse coordinates across modes and inputs across implementations.
            # Roundness derives its own weights; unity uses explicit ones.
            for mode in ("unity", "weighted", "fuzzy", "roundness"):
                routine = "pytess_roundness" if mode == "roundness" else "pytess_simple"
                args = (size, size, x, y, snr)
                if mode != "roundness":
                    args += (np.ones(count) if mode == "unity" else weights,)
                args += (0.05, "c" if mode == "fuzzy" else "s")
                functions = [getattr(fortran, routine), getattr(python, routine)]
                # Validation and warmup happen outside the timed region.
                np.testing.assert_array_equal(functions[0](*args), functions[1](*args), strict=True)
                for fn in functions:
                    fn(*args)
                times = [[], []]
                for repetition in range(options.repeat):
                    # Alternate first runner to reduce consistent ordering bias
                    # from cache state or CPU frequency changes.
                    for index in ((0, 1) if repetition % 2 == 0 else (1, 0)):
                        start = perf_counter()
                        functions[index](*args)
                        times[index].append(perf_counter() - start)
                ft, pt = (median(t) for t in times)
                # Keep every sample so readers can assess variability, not just
                # the displayed median. JSON durations are in seconds.
                rows.append(dict(shape=[size, size], generators=count, mode=mode,
                                 fortran_seconds=times[0], python_seconds=times[1],
                                 fortran_median=ft, python_median=pt, speedup=ft / pt))
                print(f"{size:4}x{size:<4} {count:10} {mode:10} {ft*1000:11.3f} {pt*1000:11.3f} {ft/pt:9.2f}x",
                      flush=True)
    # The source hash ties results to the measured implementation; environment
    # metadata helps distinguish code changes from compiler/platform effects.
    report = dict(python=sys.version, numpy=np.__version__, platform=platform.platform(),
                  processor=platform.processor(), seed=options.seed,
                  implementation_sha256=hashlib.sha256((ROOT / "bdsf" / "_tessellation.py").read_bytes()).hexdigest(),
                  compiler=subprocess.check_output(["gfortran", "--version"], text=True).splitlines()[0],
                  fortran_flags="-O3 -fno-fast-math", results=rows)
    if options.output:
        with open(options.output, "w") as stream:
            json.dump(report, stream, indent=2)
            stream.write("\n")
    if options.min_speedup is not None and any(row["speedup"] < options.min_speedup for row in rows):
        raise SystemExit("Performance target missed; see the per-case results above")


if __name__ == "__main__":
    main()

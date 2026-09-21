"""Require exact Python/Fortran equality for deterministic and random cases.

Run: python test/check_tessellation.py --random-cases 200
Requires NumPy, a C compiler and gfortran (plus Meson for NumPy 2).
"""

import argparse
from collections import Counter
from time import perf_counter

import numpy as np

from tessellation_reference import implementations


def cases(seed=20260921, random_cases=200, large_cases=4):
    """Yield (name, routine, arguments) for exact comparisons with Fortran.

    Each image layout exercises two simple weight choices, hard assignment,
    several fuzzy tolerances, and roundness. Random images are seeded; larger
    rectangular images also exercise partial blocks. Final cases target weight
    extremes and floating-point boundaries that ordinary random inputs miss.
    """
    rng = np.random.default_rng(seed)
    fixed = [
        ("one-generator", (17, 23), [8.5], [12.5]),
        ("symmetric-ties", (33, 35), [8., 26., 8., 26.], [9., 9., 27., 27.]),
        ("duplicate-generators", (17, 23), [4., 4., 15.], [6., 6., 19.]),
        ("outside-image", (13, 19), [-10., 40., 7.], [5., 5., 14.]),
        ("single-row", (1, 101), [1., 1., 1.], [5., 50., 95.]),
        ("single-column", (101, 1), [5., 50., 95.], [1., 1., 1.]),
        ("tile-boundaries", (129, 131), [1., 64., 65., 128.], [1., 65., 64., 131.]),
        ("near-ties", (65, 67), [16., np.nextafter(50., 51.)], [34., 34.]),
        ("empty-tile", (13, 19), [3., 11., 1000.], [4., 15., 1000.]),
        ("singleton-tile", (9, 9), [1., 1., 2., 7.], [1., 2., 1., 7.]),
        ("singleton-edge-block", (65, 65), [65., 65., 64., 3.], [65., 64., 65., 3.]),
        ("near-equal-weights", (67, 71), [12., 55.], [35., 35.]),
    ]
    for index in range(random_cases):
        shape = tuple(int(x) for x in rng.integers(3, 140, size=2))
        count = int(rng.integers(1, 40))
        x = rng.uniform(-2, shape[0] + 2, count)
        y = rng.uniform(-2, shape[1] + 2, count)
        fixed.append(("random-" + str(index), shape, x, y))
    for index in range(large_cases):
        shape = (512 + index * 17, 513 + index * 13)
        count = (2, 8, 32, 64)[index % 4]
        x = rng.uniform(1, shape[0], count)
        y = rng.uniform(1, shape[1], count)
        fixed.append(("large-" + str(index), shape, x, y))
    for name, (n, m), x, y in fixed:
        x, y = np.asarray(x), np.asarray(y)
        snr = np.full(x.size, 10.)
        for kind, weights in [("unity", np.ones(x.size)),
                              ("weighted", rng.uniform(0.2, 3., x.size))]:
            for code, eps in [("s", 0.05), ("c", 0.), ("c", 0.05), ("c", 1.), ("c", -0.1)]:
                yield name + "/" + kind + "/" + code + "/" + str(eps), "pytess_simple", (
                    n, m, x, y, snr, weights, eps, code)
        yield name + "/roundness", "pytess_roundness", (n, m, x, y, snr, 0.05, "s")
    for weights in ([1., 1.], [1.e-12, 1.e12], [1.e12, 1.e-12]):
        yield "weight-range-" + str(weights), "pytess_simple", (
            71, 69, [1., 71.], [1., 69.], [1., 1.], weights, 0.05, "s")
    yield "single-pixel", "pytess_simple", (1, 1, [1.], [1.], [1.], [1.], 0.05, "s")
    for weight in (np.nextafter(1., 0.), np.nextafter(1., 2.)):
        # Adjacent representable weights/tolerances probe rounding-sensitive
        # decisions without introducing an approximate comparison tolerance.
        for eps in (np.nextafter(0., 1.), np.nextafter(0.05, 0.), np.nextafter(0.05, 1.)):
            for code in ("s", "c"):
                yield "near-equal-distances", "pytess_simple", (
                    71, 69, [18., 54.], [35., 35.], [1., 1.], [1., weight], eps, code)


def check_case(python, fortran, name, routine, args):
    """Require identical values, dtype and shape, plus Fortran output layout.

    Both modules receive the same arguments. Assertion errors carry the case
    name, allowing a failing layout to be located in ``cases``.
    """
    expected = getattr(fortran, routine)(*args)
    actual = getattr(python, routine)(*args)
    np.testing.assert_array_equal(actual, expected, err_msg=name, strict=True)
    assert actual.flags.f_contiguous, name


def check_roundness_reductions(python, fortran):
    """Compare centroids and inverse radii directly with the Fortran helper.

    Duplicate generators create an empty tile. Matching its NaNs and the
    nonempty tiles' exact reductions checks behavior that final labels alone
    could conceal. The Fortran helper fills the supplied arrays in place.
    """
    n, m = 131, 137
    x, y = np.array([10., 25., 78., 78., 120.]), np.array([16., 110., 65., 65., 100.])
    labels = fortran.pytess_simple(n, m, x, y, np.ones(5), np.ones(5), 0.05, "s")
    expected_x, expected_y, radius, factor = (np.zeros(5) for _ in range(4))
    fortran.tile_roundness(labels, x, y, factor, radius, expected_x, expected_y)
    flat = labels.ravel(order="F").astype(np.intp) - 1
    with np.errstate(invalid="ignore"):
        _, cx, cy = python._centroids(flat, n, m, 5)
    np.testing.assert_array_equal(cx, expected_x)
    np.testing.assert_array_equal(cy, expected_y)
    dx = cx[flat].reshape(m, n) - np.arange(1, n + 1)[None, :]
    dy = cy[flat].reshape(m, n) - np.arange(1, m + 1)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        radii = np.sqrt(dx * dx + dy * dy).ravel()
        actual_radius = 1.0 / (np.bincount(flat, weights=radii, minlength=5)
                               / np.bincount(flat, minlength=5))
    np.testing.assert_array_equal(actual_radius, radius)


def check_invalid_inputs(python):
    """Check explicit ValueError handling and return the number of cases.

    Do not call the Fortran oracle on these inputs: some cause undefined labels
    or invalid memory accesses rather than a meaningful reference result.
    """
    simple = (4, 5, [1.], [1.], [1.], [1.], 0.05, "s")
    invalid = []
    for index, value in [(0, 0), (2, []), (3, [np.nan]), (5, [0.]),
                         (5, [-1.]), (5, [np.inf]), (6, np.nan), (7, "z")]:
        args = list(simple)
        args[index] = value
        invalid.append((python.pytess_simple, args))
    invalid.append((python.pytess_roundness, (4, 5, [1.], [1.], [1.], 0.05, "c")))
    invalid.append((python.pytess_roundness, (1, 1, [1.], [1.], [1.], 0.05, "s")))
    for fn, args in invalid:
        try:
            fn(*args)
        except ValueError:
            continue
        raise AssertionError("Expected invalid input to be rejected: " + repr(args))
    return len(invalid)


def main():
    """Run the comparison CLI, report progress, and stop at the first failure.

    Default output summarizes groups and modes. ``--verbose`` prints every
    comparison; failure context is printed in either mode before re-raising.
    Reference setup time is reported separately from time spent checking.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-cases", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260921)
    parser.add_argument("--large-cases", type=int, default=4)
    parser.add_argument("-v", "--verbose", action="store_true", help="Print every comparison")
    options = parser.parse_args()
    print("Tessellation equivalence check (exact equality, no tolerance)", flush=True)
    print(f"Seed: {options.seed}; random images: {options.random_cases}; large images: {options.large_cases}", flush=True)
    print("Loading Python implementation and building/loading the Fortran reference...", flush=True)
    started = perf_counter()
    python, fortran = implementations()
    print(f"Implementations ready ({perf_counter() - started:.2f}s)", flush=True)
    count = 0
    totals = Counter()
    phase = None
    phase_count = 0
    checked_at = perf_counter()
    for name, routine, args in cases(options.seed, options.random_cases, options.large_cases):
        current_phase = ("Randomized images" if name.startswith("random-") else
                         "Large images" if name.startswith("large-") else
                         "Deterministic and edge cases")
        if current_phase != phase:
            if phase is not None:
                print(f"  PASS: {phase_count} comparisons", flush=True)
            print(f"\n{current_phase}", flush=True)
            phase, phase_count = current_phase, 0
        detail = f"{name}: {args[0]}x{args[1]}, {len(args[2])} generators, {routine}, code={args[-1]}, eps={args[-2]}"
        if options.verbose:
            print(f"  Checking {detail}", flush=True)
        try:
            check_case(python, fortran, name, routine, args)
        except Exception:
            print(f"  FAIL: {detail}", flush=True)
            raise
        count += 1
        phase_count += 1
        mode = "Roundness" if routine == "pytess_roundness" else "Simple hard" if args[-1] == "s" else "Simple fuzzy"
        totals[mode] += 1
        if options.verbose:
            print("    PASS", flush=True)
        elif phase_count % 250 == 0:
            print(f"  {phase_count} comparisons passed; latest: {name}", flush=True)
    print(f"  PASS: {phase_count} comparisons", flush=True)
    print("\nChecking invalid inputs...", flush=True)
    invalid = check_invalid_inputs(python)
    print(f"  PASS: {invalid} invalid inputs rejected", flush=True)
    print("Checking roundness centroids and inverse mean radii...", flush=True)
    check_roundness_reductions(python, fortran)
    print("  PASS: exact intermediate results, including empty-tile NaNs", flush=True)
    print(f"\nPASS: {count} exact array comparisons in {perf_counter() - checked_at:.2f}s", flush=True)
    for mode, total in totals.items():
        print(f"  {mode}: {total}")
    print(f"  Invalid-input checks: {invalid}")
    print(f"  Total elapsed time (including reference setup): {perf_counter() - started:.2f}s")


if __name__ == "__main__":
    main()

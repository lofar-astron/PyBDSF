"""Require exact Python/Fortran equality for deterministic and random cases.

Run: python test/check_tessellation.py --random-cases 200
Requires NumPy, a C compiler and gfortran (plus Meson for NumPy 2).
"""

import argparse

import numpy as np

from tessellation_reference import implementations


def cases(seed=20260921, random_cases=200, large_cases=4):
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
        for eps in (np.nextafter(0., 1.), np.nextafter(0.05, 0.), np.nextafter(0.05, 1.)):
            for code in ("s", "c"):
                yield "near-equal-distances", "pytess_simple", (
                    71, 69, [18., 54.], [35., 35.], [1., 1.], [1., weight], eps, code)


def check_case(python, fortran, name, routine, args):
    expected = getattr(fortran, routine)(*args)
    actual = getattr(python, routine)(*args)
    np.testing.assert_array_equal(actual, expected, err_msg=name, strict=True)
    assert actual.flags.f_contiguous, name


def check_roundness_reductions(python, fortran):
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--random-cases", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260921)
    parser.add_argument("--large-cases", type=int, default=4)
    options = parser.parse_args()
    python, fortran = implementations()
    count = 0
    for name, routine, args in cases(options.seed, options.random_cases, options.large_cases):
        check_case(python, fortran, name, routine, args)
        count += 1
    invalid = check_invalid_inputs(python)
    check_roundness_reductions(python, fortran)
    print(f"PASS: {count} exact array comparisons; {invalid} invalid-input checks; exact roundness reductions")


if __name__ == "__main__":
    main()

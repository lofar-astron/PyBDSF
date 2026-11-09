"""Check that two PyBDSF output directories contain identical FITS files.

Usage: compare_fits.py <dir1> <dir2>

All FITS files (images and catalogs) found in either directory, including
subdirectories, are compared exactly: the data of every HDU must be equal
(NaNs in the same places). Headers are not compared, since they contain dates
and version numbers. Exit status 0 if all files are identical (and there is
at least one).
"""
import os
import sys

import numpy as N
from astropy.io import fits


def fits_files(top):
    found = set()
    for root, _, files in os.walk(top):
        for name in files:
            if name.endswith('.fits'):
                found.add(os.path.relpath(os.path.join(root, name), top))
    return found


def same_data(a, b):
    if a is None or b is None:
        return a is None and b is None
    if a.dtype.names:
        if a.dtype.names != b.dtype.names or a.shape != b.shape:
            return False
        return all(same_data(N.asarray(a[name]), N.asarray(b[name])) for name in a.dtype.names)
    if a.shape != b.shape:
        return False
    if a.dtype.kind in 'fc':
        return N.array_equal(a, b, equal_nan=True)
    return N.array_equal(a, b)


def main(dir1, dir2):
    files1, files2 = fits_files(dir1), fits_files(dir2)
    if not files1 and not files2:
        print(f'no FITS files found in {dir1} and {dir2}')
        return 1
    ok = True
    for name in sorted(files1 ^ files2):
        print(f'only in one directory: {name}')
        ok = False
    for name in sorted(files1 & files2):
        with fits.open(os.path.join(dir1, name)) as h1, fits.open(os.path.join(dir2, name)) as h2:
            if len(h1) != len(h2):
                print(f'{name}: different number of HDUs')
                ok = False
                continue
            for i, (a, b) in enumerate(zip(h1, h2)):
                if not same_data(a.data, b.data):
                    print(f'{name}: HDU {i} differs')
                    ok = False
    print(f'{len(files1 & files2)} FITS files compared: '
          f'{"all identical" if ok else "DIFFERENCES FOUND"}')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:3]))

#!/bin/bash
# Compare the pure-Python natgrid (working tree) with the C extension
# natgridmodule of a baseline revision (default: master).
#
# Usage: run_comparison.sh <work dir> [baseline git revision] [repeats]
#
# Three versions are compared:
#   orig : the baseline
#   fix  : the baseline, with contiguous coordinate arrays passed to
#          natgridmodule. In the baseline, psf_vary passes the strided views
#          from numpy.transpose(psfcoords), and natgridmodule (which ignores
#          strides) reads interleaved x and y coordinates. `fix` computes what
#          was intended, and is the reference for the outputs of `new`.
#   new  : the working tree, with bdsf/nat (natgridmodule) removed
#
# Steps:
#   1. Build overlay packages pkg_orig, pkg_fix and pkg_new, using
#      ../rmsimage-optimization/setup_overlays.sh
#   2. Run test/tnatgrid.py with the new version
#   3. Run test/tbdsf_process_image.py with all versions, and compare the FITS
#      outputs of new and fix
#   4. Time benchmark/run_pybdsf.py with all versions (CPU time and peak
#      memory from /usr/bin/time, wall time of each step from print_timing),
#      and compare the FITS outputs of new and fix
set -e
WORK=${1:?usage: $0 <work dir> [baseline git revision] [repeats]}
REV=${2:-master}
REPS=${3:-3}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel)
PY=$REPO/venv/bin/python
FITS=$REPO/benchmark/NGC6048_9-MFS-image.fits
VERSIONS="orig fix new"

"$HERE/../rmsimage-optimization/setup_overlays.sh" "$WORK" "$REV"
rm -rf "$WORK/pkg_new/bdsf/nat"
rm -rf "$WORK/pkg_fix"
cp -r "$WORK/pkg_orig" "$WORK/pkg_fix"
sed -i 's/x, y = N.transpose(psfcoords)/x, y = N.ascontiguousarray(N.transpose(psfcoords))/' \
    "$WORK/pkg_fix/bdsf/psf_vary.py"
[ "$(grep -c 'N.ascontiguousarray(N.transpose(psfcoords))' "$WORK/pkg_fix/bdsf/psf_vary.py")" = 3 ]

compare() {  # compare <dir ref> <dir new>: max relative difference of FITS data
    $PY -I -c "
import pathlib
import numpy as np
from astropy.io import fits
ref, new = pathlib.Path('$1'), pathlib.Path('$2')
files = sorted(p.relative_to(ref) for p in ref.rglob('*') if p.suffix.lower() == '.fits')
identical, worst = 0, []
for f in files:
    if not (new / f).exists():
        worst.append((float('inf'), str(f), 'missing')); continue
    same, maxrel = True, 0.0
    with fits.open(ref / f) as a, fits.open(new / f) as b:
        for ha, hb in zip(a, b):
            da, db = ha.data, hb.data
            if da is None: continue
            if isinstance(ha, fits.BinTableHDU):
                for name in da.names:
                    ca, cb = np.asarray(da[name]), np.asarray(db[name])
                    if ca.dtype.kind in 'fc':
                        same &= np.array_equal(ca, cb, equal_nan=True)
                        d = np.abs(ca - cb) / np.maximum(np.abs(ca), 1e-30)
                        maxrel = max(maxrel, float(np.nanmax(d, initial=0)))
                    else:
                        same &= np.array_equal(ca, cb)
            else:
                same &= np.array_equal(da, db, equal_nan=True)
                scale = np.nanmax(np.abs(da), initial=0) or 1.0
                maxrel = max(maxrel, float(np.nanmax(np.abs(da - db), initial=0)) / scale)
    identical += bool(same)
    if not same: worst.append((maxrel, str(f), ''))
print('    %d FITS files compared, %d identical' % (len(files), identical))
for m, f, s in sorted(worst, reverse=True)[:8]:
    print('      %-60s %s' % (f, s or 'max rel. diff %.2e' % m))"
}

echo "=== unit test (new)"
(cd "$REPO/test" && PYTHONPATH="$WORK/pkg_new" $PY tnatgrid.py > "$WORK/tnatgrid.log" 2>&1) \
    && echo "tnatgrid.py: passed" || echo "tnatgrid.py: FAILED (see $WORK/tnatgrid.log)"

echo "=== test/tbdsf_process_image.py"
for v in $VERSIONS; do
    rm -rf "$WORK/test_$v"
    mkdir "$WORK/test_$v"
    cp -r "$REPO"/test/{tbdsf_process_image.py,compare_outputs.py,tbdsf_process_image.in,tbdsf_process_image.in_fits,reference_outputs} "$WORK/test_$v/"
    (cd "$WORK/test_$v" && PYTHONPATH="$WORK/pkg_$v" $PY tbdsf_process_image.py > run.log 2>&1) \
        && status=passed || status="FAILED (see $WORK/test_$v/run.log)"
    echo "$v: $status"
done
echo "  new vs fix:"
compare "$WORK/test_fix/tbdsf_process_image.in_fits_pybdsf" "$WORK/test_new/tbdsf_process_image.in_fits_pybdsf"

echo "=== run_pybdsf.py"
for i in $(seq "$REPS"); do
    for v in $VERSIONS; do
        rm -rf "$WORK/run_$v"
        mkdir "$WORK/run_$v"
        ln -s "$FITS" "$WORK/run_$v/"
        sed 's/clobber = True)/clobber = True, print_timing = True)/' \
            "$REPO/benchmark/run_pybdsf.py" > "$WORK/run_$v/run_pybdsf.py"
        (cd "$WORK/run_$v" && PYTHONPATH="$WORK/pkg_$v" /usr/bin/time -v $PY run_pybdsf.py > out.log 2> time.log)
        cpu=$(awk -F': ' '/User time/ {u=$2} /System time/ {s=$2} END {printf "%.1f", u + s}' "$WORK/run_$v/time.log")
        rss=$(awk -F': ' '/Maximum resident/ {printf "%.0f", $2 / 1024}' "$WORK/run_$v/time.log")
        echo "$v run $i: CPU ${cpu} s, peak RSS ${rss} MB, $(grep -E '^ +(Op_psf_vary|Total) :' "$WORK/run_$v/out.log" | tr -s ' ' | tr '\n' ' ')"
    done
done
echo "  new vs fix:"
compare "$WORK/run_fix/NGC6048_9-MFS-image_pybdsf" "$WORK/run_new/NGC6048_9-MFS-image_pybdsf"
echo "  new vs orig:"
compare "$WORK/run_orig/NGC6048_9-MFS-image_pybdsf" "$WORK/run_new/NGC6048_9-MFS-image_pybdsf"

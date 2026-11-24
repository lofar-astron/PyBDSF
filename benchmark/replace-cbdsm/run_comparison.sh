#!/bin/bash
# Compare the pure-Python MGFunction (working tree) with the C++ `_cbdsm` of
# a baseline revision (default: master).
#
# Usage: run_comparison.sh <work dir> [baseline git revision] [repeats]
#
# Steps:
#   1. Build overlay packages: pkg_orig (baseline) and pkg_new (working tree,
#      with `_cbdsm` removed), using ../rmsimage-optimization/setup_overlays.sh
#   2. Run test/tmgfunction.py and test/tmulti_proc.py with the new version
#   3. Run test/tbdsf_process_image.py with both versions (timed), and check
#      that all their FITS outputs (images and catalogs) are identical
#   4. Time benchmark/run_pybdsf.py with both versions (with print_timing),
#      and check that their outputs are identical
set -e
WORK=${1:?usage: $0 <work dir> [baseline git revision] [repeats]}
REV=${2:-master}
REPS=${3:-3}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel)
PY=$REPO/venv/bin/python
FITS=$REPO/benchmark/NGC6048_9-MFS-image.fits

"$HERE/../rmsimage-optimization/setup_overlays.sh" "$WORK" "$REV"
rm -f "$WORK"/pkg_new/bdsf/_cbdsm*.so

compare() {  # compare <dir orig> <dir new>: all FITS files must be identical
    $PY -I -c "
import pathlib, sys
from astropy.io.fits import FITSDiff
orig, new = pathlib.Path('$1'), pathlib.Path('$2')
files = sorted(p.relative_to(orig) for p in orig.rglob('*') if p.suffix.lower() == '.fits')
differ = [str(f) for f in files
          if not (new / f).exists()
          or not FITSDiff(str(orig / f), str(new / f), rtol=0, ignore_keywords=['*']).identical]
print('    %d FITS files compared, %d differ %s' % (len(files), len(differ), differ or ''))"
}

echo "=== unit tests (new)"
for t in tmgfunction.py tmulti_proc.py; do
    (cd "$REPO/test" && PYTHONPATH="$WORK/pkg_new" $PY "$t" > "$WORK/$t.log" 2>&1) \
        && echo "$t: passed" || echo "$t: FAILED (see $WORK/$t.log)"
done

echo "=== test/tbdsf_process_image.py"
for i in $(seq "$REPS"); do
    for v in orig new; do
        rm -rf "$WORK/test_$v"
        mkdir "$WORK/test_$v"
        cp -r "$REPO"/test/{tbdsf_process_image.py,compare_outputs.py,tbdsf_process_image.in,tbdsf_process_image.in_fits,reference_outputs} "$WORK/test_$v/"
        start=$(date +%s.%N)
        (cd "$WORK/test_$v" && PYTHONPATH="$WORK/pkg_$v" $PY tbdsf_process_image.py > run.log 2>&1) \
            && status=passed || status="FAILED (see $WORK/test_$v/run.log)"
        echo "$v run $i: $status, $(echo "$(date +%s.%N) - $start" | bc) s"
    done
done
compare "$WORK/test_orig/tbdsf_process_image.in_fits_pybdsf" "$WORK/test_new/tbdsf_process_image.in_fits_pybdsf"

echo "=== run_pybdsf.py"
for i in $(seq "$REPS"); do
    for v in orig new; do
        rm -rf "$WORK/run_$v"
        mkdir "$WORK/run_$v"
        ln -s "$FITS" "$WORK/run_$v/"
        sed 's/clobber = True)/clobber = True, print_timing = True)/' \
            "$REPO/benchmark/run_pybdsf.py" > "$WORK/run_$v/run_pybdsf.py"
        (cd "$WORK/run_$v" && PYTHONPATH="$WORK/pkg_$v" $PY run_pybdsf.py > out.log 2>&1)
        echo "$v run $i: $(grep -E '^ +(Op_gausfit|Op_wavelet_atrous|Total) :' "$WORK/run_$v/out.log" | tr -s ' ' | tr '\n' ' ')"
    done
done
compare "$WORK/run_orig/NGC6048_9-MFS-image_pybdsf" "$WORK/run_new/NGC6048_9-MFS-image_pybdsf"

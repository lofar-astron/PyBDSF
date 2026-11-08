#!/bin/bash
# Compare the Python tessellation (working tree) with the Fortran extension
# `_pytesselate` (baseline revision, default: master).
#
# Usage: run_comparison.sh <work dir> [baseline git revision]
#
# Steps:
#   1. Build overlay packages: pkg_orig (baseline) and pkg_new (working tree),
#      using ../rmsimage-optimization/setup_overlays.sh. The installed
#      package in venv must still contain `_pytesselate`.
#   2. Run test/ttessellate.py with the new version
#   3. Run test/tbdsf_process_image.py with both versions, and check that all
#      FITS outputs are identical
#   4. Run run_psf_vary.py on NGC6048 with both versions, for the simple
#      (psf_itess_method=0) and the roundness (1) tessellation, and check
#      that all FITS outputs are identical
set -e
WORK=${1:?usage: $0 <work dir> [baseline git revision]}
REV=${2:-master}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel)
PY=$REPO/venv/bin/python
FITS=$REPO/benchmark/NGC6048_9-MFS-image.fits

"$HERE/../rmsimage-optimization/setup_overlays.sh" "$WORK" "$REV"

echo "=== test/ttessellate.py (new)"
(cd "$REPO/test" && PYTHONPATH="$WORK/pkg_new" $PY ttessellate.py)

echo "=== test/tbdsf_process_image.py"
for v in orig new; do
    rm -rf "$WORK/test_$v"
    mkdir "$WORK/test_$v"
    cp -r "$REPO"/test/{tbdsf_process_image.py,compare_outputs.py,tbdsf_process_image.in,tbdsf_process_image.in_fits,reference_outputs} "$WORK/test_$v/"
    (cd "$WORK/test_$v" && PYTHONPATH="$WORK/pkg_$v" $PY tbdsf_process_image.py > run.log 2>&1 \
        && echo "$v: test passed" || echo "$v: test FAILED (see $WORK/test_$v/run.log)")
done
PYTHONPATH="$WORK/pkg_new" $PY "$HERE/compare_fits.py" \
    "$WORK/test_orig/tbdsf_process_image.in_fits_pybdsf" "$WORK/test_new/tbdsf_process_image.in_fits_pybdsf"

echo "=== run_psf_vary.py (NGC6048)"
for method in 0 1; do
    for v in orig new; do
        D="$WORK/run_${v}_m$method"
        rm -rf "$D"
        mkdir "$D"
        ln -s "$FITS" "$D/"
        (cd "$D" && PYTHONPATH="$WORK/pkg_$v" $PY "$HERE/run_psf_vary.py" "$(basename "$FITS")" $method > out.log 2>&1; \
            echo "method $method, $v: $(grep 'running time' out.log)")
    done
    PYTHONPATH="$WORK/pkg_new" $PY "$HERE/compare_fits.py" "$WORK/run_orig_m$method" "$WORK/run_new_m$method"
done

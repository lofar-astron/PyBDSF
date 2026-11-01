#!/bin/bash
# Run all benchmarks and comparisons used for the rmsimage/bstat optimization.
#
# Usage: run_comparison.sh <work dir> [git revision of the baseline, default HEAD]
#
# Steps:
#   1. Build overlay packages (baseline revision vs working tree)
#   2. bstat speed/accuracy comparison
#   3. rms_mean_map timings + map comparison (normal and 2x2 tiled image)
#   4. End-to-end test (test/tbdsf_process_image.py) with both versions, and a
#      direct comparison of their outputs at rtol 1e-3 (catalogs included)
#   5. Full benchmark/run_pybdsf.py with both versions
set -e
WORK=${1:?usage: $0 <work dir> [git revision]}
REV=${2:-HEAD}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel)
PY=$REPO/venv/bin/python
FITS=$REPO/benchmark/NGC6048_9-MFS-image.fits

"$HERE/setup_overlays.sh" "$WORK" "$REV"
cd "$WORK"

echo "=== bstat"
$PY -I "$HERE/compare_bstat.py" "$WORK/pkg_orig" "$WORK/pkg_new" "$FITS"

echo "=== rms_mean_map"
for v in orig new; do
    $PY -I "$HERE/bench_rms_mean_map.py" "$WORK/pkg_$v" $v "$FITS" "$WORK/maps_$v.npz"
    $PY -I "$HERE/bench_rms_mean_map.py" "$WORK/pkg_$v" $v "$FITS" "$WORK/maps_tiled_$v.npz" --tiled
done
$PY "$HERE/compare_maps.py" "$WORK/maps_orig.npz" "$WORK/maps_new.npz"
$PY "$HERE/compare_maps.py" "$WORK/maps_tiled_orig.npz" "$WORK/maps_tiled_new.npz"

echo "=== end-to-end test"
for v in orig new; do
    rm -rf "$WORK/test_$v"
    mkdir "$WORK/test_$v"
    cp -r "$REPO"/test/{tbdsf_process_image.py,compare_outputs.py,tbdsf_process_image.in,tbdsf_process_image.in_fits,reference_outputs} "$WORK/test_$v/"
    (cd "$WORK/test_$v" && PYTHONPATH="$WORK/pkg_$v" $PY tbdsf_process_image.py > run.log 2>&1 \
        && echo "$v: test passed" || echo "$v: test FAILED (see $WORK/test_$v/run.log)")
done
(cd "$WORK/test_new" && PYTHONPATH="$WORK/pkg_new" $PY -c "
import filecmp
from compare_outputs import compare_results
ok = compare_results(filecmp.dircmp('$WORK/test_orig/tbdsf_process_image.in_fits_pybdsf',
                                    'tbdsf_process_image.in_fits_pybdsf'),
                     0.001, 1, check_images=True, check_catalogs=True, check_other=False)
print('orig vs new outputs agree at rtol 1e-3:', ok)")

echo "=== run_pybdsf.py"
for v in orig new; do
    rm -rf "$WORK/run_$v"
    mkdir "$WORK/run_$v"
    ln -s "$FITS" "$WORK/run_$v/"
    cp "$REPO/benchmark/run_pybdsf.py" "$WORK/run_$v/"
    (cd "$WORK/run_$v" && PYTHONPATH="$WORK/pkg_$v" $PY run_pybdsf.py > out.log 2>&1; \
        echo "$v: $(grep 'running time' out.log)")
done

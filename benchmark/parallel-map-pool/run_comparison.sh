#!/bin/bash
# Compare the parallel_map rewrite (working tree) with a baseline revision
# (default: master).
#
# Usage: run_comparison.sh <work dir> [baseline git revision]
#
# Steps:
#   1. Build overlay packages: pkg_orig (baseline) and pkg_new (working tree),
#      using setup_overlays.sh
#   2. Run test/tmulti_proc.py with the new version
#   3. Run test/tbdsf_process_image.py with both versions, and compare their
#      outputs at rtol 1e-3 (catalogs included)
#   4. Time benchmark/run_pybdsf.py with both versions and print the key
#      results from the logs
set -e
WORK=${1:?usage: $0 <work dir> [baseline git revision]}
REV=${2:-master}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(git -C "$HERE" rev-parse --show-toplevel)
PY=$REPO/venv/bin/python
FITS=$REPO/benchmark/NGC6048_9-MFS-image.fits

"$HERE/setup_overlays.sh" "$WORK" "$REV"

echo "=== test/tmulti_proc.py (new)"
(cd "$REPO/test" && PYTHONPATH="$WORK/pkg_new" $PY -W error::DeprecationWarning tmulti_proc.py)

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
        echo "$v: $(grep 'running time' out.log)"; \
        grep -h -E "Min/max values of background rms|Number of islands found|Total flux density in model" \
            NGC6048_9-MFS-image.fits.pybdsf.log | head -3 | sed 's/^.*USERINFO: /    /')
done

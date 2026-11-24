#!/bin/bash
# Build two copies of the installed bdsf package (which contains the compiled
# extensions), with the pure-Python sources replaced by:
#   pkg_orig/bdsf : the sources at a git revision (default: HEAD)
#   pkg_new/bdsf  : the sources in the working tree
# The benchmark and comparison scripts take one of these directories as the
# first entry on sys.path, so both versions can be compared side by side
# without rebuilding PyBDSF.
#
# Usage: setup_overlays.sh <output dir> [git revision]
set -e
OUT=${1:?usage: $0 <output dir> [git revision]}
REV=${2:-HEAD}
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
INSTALLED=$("$REPO/venv/bin/python" -I -c 'import bdsf, os; print(os.path.dirname(bdsf.__file__))' 2>/dev/null)

mkdir -p "$OUT"
for v in orig new; do
    rm -rf "$OUT/pkg_$v"
    mkdir -p "$OUT/pkg_$v"
    cp -r "$INSTALLED" "$OUT/pkg_$v/"
    # Shared libraries bundled with the wheel (e.g. libgfortran)
    [ -d "$INSTALLED.libs" ] && cp -r "$INSTALLED.libs" "$OUT/pkg_$v/"
done
rm -rf "$OUT/src_orig"
mkdir -p "$OUT/src_orig"
git -C "$REPO" archive "$REV" bdsf | tar -x -C "$OUT/src_orig"
cp "$OUT"/src_orig/bdsf/*.py "$OUT/pkg_orig/bdsf/"
cp "$REPO"/bdsf/*.py "$OUT/pkg_new/bdsf/"
echo "Created $OUT/pkg_orig ($REV) and $OUT/pkg_new (working tree)"

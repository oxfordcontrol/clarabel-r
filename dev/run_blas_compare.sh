#!/bin/bash
# Run the SDP benchmark under the reference BLAS and then under vecLib,
# restoring the reference BLAS on any exit path (success, error, or signal).
#
# The symlink is system-wide for newly started R processes, so this must not be
# left pointing at vecLib.
set -u

LIB=/Library/Frameworks/R.framework/Versions/4.6/Resources/lib
SCRATCH=/private/tmp/claude-503/-Users-naras-GitHub-clarabel-r/826c161b-9109-489f-b745-15dabc1b237a/scratchpad
ORIG=$(readlink "$LIB/libRblas.dylib")

restore() {
  ln -sf libRblas.0.dylib "$LIB/libRblas.dylib"
  echo "=== restored: libRblas.dylib -> $(readlink "$LIB/libRblas.dylib") (was $ORIG at start) ==="
}
trap restore EXIT INT TERM

echo "=== reference BLAS run ==="
ln -sf libRblas.0.dylib "$LIB/libRblas.dylib"
Rscript "$SCRATCH/blas_compare.R" "$SCRATCH/res_reference.rds"

echo
echo "=== vecLib run ==="
ln -sf libRblas.vecLib.dylib "$LIB/libRblas.dylib"
Rscript -e 'cat("verify BLAS in use:", extSoftVersion()[["BLAS"]], "\n")'
Rscript "$SCRATCH/blas_compare.R" "$SCRATCH/res_veclib.rds"

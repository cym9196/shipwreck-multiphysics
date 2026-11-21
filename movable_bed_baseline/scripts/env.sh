#!/usr/bin/env bash
# Source before invoking OpenFOAM; outputs are isolated in this project.
project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
if [[ -n "${OPENFOAM_BASHRC:-}" ]]; then
    source "$OPENFOAM_BASHRC" || return $?
fi
if ! command -v wmake >/dev/null 2>&1; then
    echo 'Load the OpenFOAM environment or set OPENFOAM_BASHRC to its etc/bashrc.' >&2
    return 1
fi
export SHIPWRECK_BASE="$project_root"
export FOAM_USER_APPBIN="$project_root/build/bin"
export FOAM_USER_LIBBIN="$project_root/build/lib"
export FOAM_MODULE_APPBIN="$FOAM_USER_APPBIN"
export FOAM_MODULE_LIBBIN="$FOAM_USER_LIBBIN"
export PATH="$FOAM_USER_APPBIN:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OF_VERSION=$(wmake -version | cut -c1-3)


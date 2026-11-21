#!/usr/bin/env bash
scenario_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
baseline=$(realpath "${1:?Pass the movable-bed baseline directory}")
# The installed OpenFOAM startup script expects unset shell variables.
case $- in *u*) scenario_nounset=true; set +u;; *) scenario_nounset=false;; esac
scenario_env_status=0
source "$baseline/scripts/env.sh" || scenario_env_status=$?
if $scenario_nounset; then set -u; fi
if (( scenario_env_status != 0 )); then return "$scenario_env_status"; fi
export SHIPWRECK_BASE="$baseline"
export FOAM_USER_APPBIN="$scenario_root/build/bin"
export FOAM_MODULE_APPBIN="$FOAM_USER_APPBIN"
export PATH="$FOAM_USER_APPBIN:$PATH"

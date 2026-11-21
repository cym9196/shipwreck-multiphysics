#!/usr/bin/env bash
set -euo pipefail
project=$(cd -- "$(dirname -- "$0")" && pwd)
base=${1:?Pass the movable-bed baseline directory}
case_dir=$(realpath "${2:?Pass a newly prepared case}")
ranks=${3:-6}
mode=${4:-}
if [[ "$ranks" != 1 && "$ranks" != 6 ]]; then echo 'Use 1 or 6 ranks' >&2; exit 2; fi
if [[ -n "$mode" && "$mode" != --check-saved ]]; then echo 'Unknown mode' >&2; exit 2; fi
source "$project/env.sh" "$base"
cd -- "$case_dir"
test -f scenario.json
exec 9>>.run.lock
if ! flock -n 9; then
    echo "Could not lock case at $case_dir; another run or check may be active." >&2
    exit 2
fi
if [[ "$mode" != --check-saved ]]; then
    if [[ -e log.solver || -e constant/polyMesh || -e processor0 ]]; then
        echo 'Existing run detected. Prepare an independent case or diagnose the saved state.' >&2
        exit 2
    fi
    surfaceCheck -checkSelfIntersection constant/triSurface/shipwreck.stl > log.surfaceCheck 2>&1
    if ! grep -q 'Surface is closed' log.surfaceCheck || ! grep -q 'Surface is not self-intersecting' log.surfaceCheck; then
        echo 'Hull surface check failed. See log.surfaceCheck and any selfInterPoints.obj.' >&2
        exit 3
    fi
    blockMesh > log.blockMesh 2>&1
    checkMesh > log.checkMesh 2>&1
    grep -q 'Mesh OK\.' log.checkMesh
    makeFaMesh > log.makeFaMesh 2>&1
    if [[ "$ranks" == 6 ]]; then
        decomposePar > log.decomposePar 2>&1
        mpirun --bind-to core -np 6 shipwreckSedExnerFoam -parallel > log.solver 2>&1
        reconstructPar -latestTime > log.reconstructPar 2>&1
    else
        shipwreckSedExnerFoam > log.solver 2>&1
    fi
fi
grep -q '^End$' log.solver
grep -q 'Mesh OK\.' log.checkMesh
checkMesh -latestTime -allTopology -allGeometry > log.checkMesh.final 2>&1
if ! grep -q 'Mesh OK\.' log.checkMesh.final; then
    meshAnisotropy -latestTime > log.meshAnisotropy.final 2>&1
fi
python3 "$project/check_run.py" "$case_dir"

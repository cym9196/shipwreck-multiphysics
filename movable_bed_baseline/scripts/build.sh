#!/usr/bin/env bash
source "$(dirname -- "$0")/env.sh" || exit $?
set -e -o pipefail
mkdir -p "$project_root/build/bin" "$project_root/build/lib" "$project_root/logs"
export WM_NCOMPPROCS=4
for library in projectedFiniteArea settlingModels bedloadModels criticalShieldsModels sedimentBed; do
    wmake -j4 libso "$project_root/src/sedExnerFoam/src/$library"
done
wmake -j4 libso "$project_root/src/wallFunctions/roughWallFunctions"
wmake -j4 "$project_root/src/sedExnerFoam/solver"
wmake -j2 "$project_root/src/mesh-diagnostics"
wmake -j2 libso "$project_root/src/reference-column-motion"
wmake -j2 "$project_root/src/motion-replay"
ldd "$project_root/build/bin/sedExnerFoam" > "$project_root/logs/linked_libraries.txt"
if grep -q 'not found' "$project_root/logs/linked_libraries.txt"; then
    echo 'Unresolved shared-library dependency' >&2
    exit 1
fi
python3 - "$project_root" <<'PY'
from pathlib import Path
import hashlib,json,sys
r=Path(sys.argv[1])
files=[p for p in (r/'build').glob('*/*') if p.is_file()]
record={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(r/'logs/build_hashes.json').write_text(json.dumps(record,indent=2))
print('BUILD_COMPLETE',len(files),'binary/library artifacts')
PY

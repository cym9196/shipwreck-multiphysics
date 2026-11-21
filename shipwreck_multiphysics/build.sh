#!/usr/bin/env bash
# Build the independent solver from the checkpoint-capable OpenFOAM baseline.
set -euo pipefail
project=$(cd -- "$(dirname -- "$0")" && pwd)
base=$(realpath "${1:?Pass the movable-bed baseline directory}")
mkdir -p "$project/build/bin" "$project/build/lib"
exec 9>"$project/build/.build.lock"
if ! flock -n 9; then
    echo 'Another build is using this directory.' >&2
    exit 1
fi

rust_compiler=${RUSTC:-rustc}
if [[ -z ${RUSTC:-} && -x "$project/toolchain/root/usr/bin/rustc" ]]; then
    rust_compiler="$project/toolchain/root/usr/bin/rustc"
    export LD_LIBRARY_PATH="$project/toolchain/root/usr/lib/x86_64-linux-gnu${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi
source "$project/env.sh" "$base"
python3 "$project/install_solver.py" "$base" "$project"

tools_stage=$(mktemp -d "$project/build/.tools.XXXXXX")
cleanup() {
    rm -f -- "$tools_stage/test-impact" "$tools_stage/shipwreck-case" "$tools_stage/libshipwreck_impact.a"
    rmdir -- "$tools_stage"
}
trap cleanup EXIT
"$rust_compiler" --version
"$rust_compiler" --edition 2021 -D warnings --test "$project/rust/impact.rs" -o "$tools_stage/test-impact"
"$tools_stage/test-impact"
"$rust_compiler" --edition 2021 -D warnings -C opt-level=3 -C panic=abort \
    --crate-type staticlib --crate-name shipwreck_impact "$project/rust/impact.rs" \
    -o "$tools_stage/libshipwreck_impact.a"
g++ -std=c++17 -Wall -Wextra -Wpedantic -O2 "$project/src/prepare.cpp" -o "$tools_stage/shipwreck-case"

export SHIPWRECK_IMPACT_LIB="$tools_stage"
if ! cmp -s "$tools_stage/libshipwreck_impact.a" "$project/build/lib/libshipwreck_impact.a"; then
    # wmake may not track static archives passed through EXE_LIBS.
    touch "$project/build/solver/sedExnerFoam.C"
fi
wmake -j 2 "$project/build/solver"
mv -f -- "$tools_stage/libshipwreck_impact.a" "$project/build/lib/"
mv -f -- "$tools_stage/shipwreck-case" "$tools_stage/test-impact" "$project/build/bin/"

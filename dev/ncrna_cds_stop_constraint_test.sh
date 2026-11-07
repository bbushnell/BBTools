#!/usr/bin/env bash
# Finite endpoint-search invariants; no genome input or resource construction.
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$root/calcmem.sh"
z=-Xmx512m
set +eu
setEnvironment
parseXmx "$@"
set -eu
vectorFlags=()
if [[ -n "$SIMD" ]]; then read -r -a vectorFlags <<< "$SIMD"; fi
exec java -ea -XX:+ExitOnOutOfMemoryError "$z" "${vectorFlags[@]}" -cp "$root/current" test.prok.NcrnaCdsStopConstraintTest "$@"

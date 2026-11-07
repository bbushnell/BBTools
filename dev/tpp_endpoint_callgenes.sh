#!/usr/bin/env bash
# Experimental caller; uses the ordinary BBTools environment and JVM setup.
set -euo pipefail
pushd . > /dev/null
script=${BASH_SOURCE[0]}
while [[ -h "$script" ]]; do
	cd "$(dirname "$script")"
	target=$(readlink "$(basename "$script")")
	script=$target
done
cd "$(dirname "$script")/.."
root=$(pwd)
popd > /dev/null
if [[ $# == 0 || ${1:-} == --help || ${1:-} == -h ]]; then
	cat <<'USAGE'
Experimental TPP endpoint caller. Ordinary CallGenes arguments are forwarded.
tppcdsstop=f       Enable with t to limit TPP stop searches using pre-path CDS starts.
constraintlog=FILE Optional fresh-only TSV of candidate sweeps, bounds and exclusions.
Omitted/false experiment flag preserves the ordinary caller path. No default change.
USAGE
	exit 0
fi
source "$root/calcmem.sh"
z=-Xmx8g
set +eu
setEnvironment
parseXmx "$@"
set -eu
vectorFlags=()
if [[ -n "$SIMD" ]]; then read -r -a vectorFlags <<< "$SIMD"; fi
exec java -ea -XX:+ExitOnOutOfMemoryError "$z" "${vectorFlags[@]}" -cp "$root/current" test.prok.TppEndpointCallGenes "$@"

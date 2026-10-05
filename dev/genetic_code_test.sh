#!/bin/bash
# Synthetic GeneticCode checks. Compile the full dna package first.
# Optional BBTOOLS_TEST_CLASSES names an absolute private compiled-class directory.
set -eo pipefail

script_path="${BASH_SOURCE[0]}"
while [[ -L "$script_path" ]]; do
	script_dir="$(cd -- "$(dirname -- "$script_path")" && pwd)"
	script_path="$(readlink -- "$script_path")"
	[[ "$script_path" = /* ]] || script_path="$script_dir/$script_path"
done
DIR="$(cd -- "$(dirname -- "$script_path")/.." && pwd)/"
CP="${DIR}current/"
if [[ -n "${BBTOOLS_TEST_CLASSES:-}" ]]; then
	[[ "$BBTOOLS_TEST_CLASSES" = /* && -d "$BBTOOLS_TEST_CLASSES" ]] || { echo 'BBTOOLS_TEST_CLASSES must name an existing absolute directory' >&2; exit 1; }
	CP="$BBTOOLS_TEST_CLASSES:$CP"
fi
# Canonical setup currently reads optional unset variables, so enable nounset afterwards.
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
set -u
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" dna.GeneticCodeTest

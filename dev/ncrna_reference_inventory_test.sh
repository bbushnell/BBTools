#!/bin/bash
# Native reference-inventory fixtures; no genome payloads or scientific runs.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: ncrna_reference_inventory_test.sh out=/absolute/fresh_fixture_directory
Uses fixed synthetic references, including intentionally nonexistent FASTAs.
Checks union/member positive-input counts, retained loci, source strata,
cobalamin reference-view parity and malformed metadata rejection.
Run only after compiling the complete prok and test.prok packages on Dori.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaReferenceInventoryTest "$@"

#!/bin/bash
# Native asymmetric minus-strand capture and manifest fixtures.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_window_capture_driver_test.sh out=/absolute/fresh_fixture_directory -Xmx2g -eoom
Uses the current RF00174 seed resource on two tiny synthetic FASTAs.
Checks native minus-strand coordinates, borrowed-array identity, captured bytes,
zero-window input accounting, duplicate FASTA rejection and summary semantics.
The index deliberately rejects after capture; this is not an accuracy benchmark.
Compile complete affected packages on Dori before running.
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
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.CaptureNcrnaWindowsTest "$@"

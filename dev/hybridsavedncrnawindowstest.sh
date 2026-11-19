#!/bin/bash
# Tiny saved-cache native/selected/hybrid integration and rejection fixtures.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/hybridsavedncrnawindowstest.sh out=fresh_fixture_directory t=1 -Xmx512m -eoom
Checks native model-count/full-header binding, native/canonical event tuples,
selected-method locus minima, original ambiguity and unrankable dispositions,
UTF-8 panel names, one/two-worker equality and failures before publication.
Uses literal cache/rank witnesses; no real genome, verifier or performance run.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.HybridSavedNcrnaWindowsTest "$@"

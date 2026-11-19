#!/bin/bash
# Finite member-vote and default-off candidate-window observation fixtures.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_ddl_test.sh -Xmx512m -eoom
Checks member plurality, ties, empty buckets, 2,048-bucket key storage and
buffer ownership, then observation before an intentionally rejecting index.
No genome payloads, ranking experiment or additional model alignments.
Compile the complete prok and test.prok packages on Dori before running.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	RESOLVED_SCRIPT="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$RESOLVED_SCRIPT"
done
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaDdlProfileTest "$@"
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaWindowCaptureTest "$@"

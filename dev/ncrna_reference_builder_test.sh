#!/bin/bash
set -eo pipefail
usage(){
	cat <<'HELP'
Usage: dev/ncrna_reference_builder_test.sh out=new_fixture_directory

Runs native end-to-end reference-evidence fixtures: one/two-worker row equality,
model/source accounting, unassessable intervals, opposite-strand links, zero-hit
models, source reuse, bounded samples and worker-failure propagation.
Also checks saved-table quantiles, missing/empty/mixed GFF claims, Unicode paths,
population conservation and malformed-input rejection without readable payloads.
Geometry fixtures cover transitive groups, alternative full hits, partial/outer
bridges, separate input/strand/target identities, invalid intervals and zero hits.
No scientific model searches or biological truth admission occur.
HELP
}
if [ $# -eq 0 ] || [ "$1" = -h ] || [ "$1" = --help ]; then
	usage
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	RESOLVED_SCRIPT="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$RESOLVED_SCRIPT"
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.BuildNcrnaReferenceTest "$@"

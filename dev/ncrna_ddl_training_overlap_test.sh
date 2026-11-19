#!/bin/bash
# Literal alias, nonvoting-record and unresolved-binding metadata fixtures.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_ddl_training_overlap_test.sh out=fresh_directory t=1 -Xmx512m -eoom
Checks sample-path versus distinct-assembly counts, literal UTF-8 path aliases,
kept votes versus orphans/small clusters, unresolved donors and assembly pairs,
CLI parity, complete model counts and loud malformed-input rejection.
Metadata only; all sequence paths intentionally do not exist.
Compile the complete canonical current/ tree on Dori before running.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaDdlTrainingOverlapTest "$@"

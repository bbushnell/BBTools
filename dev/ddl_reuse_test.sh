#!/bin/bash
# Finite native sketch-clear and comparison-scratch fixtures.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ddl_reuse_test.sh -Xmx512m -eoom
Checks independent sequence reuse after a raised adaptive floor, optional
count/GC/kmer arrays, retained array ownership and literal comparison counts.
No external inputs or scientific ranking run. This is not a speed benchmark.
Compile complete cardinality and test.cardinality packages on Dori first.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.cardinality.DynamicDemiLogReuseTest "$@"

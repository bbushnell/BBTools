#!/bin/bash
# Literal complete-contig equality, conservative nonmatches and malformed inputs.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_construction_assembly_compare_test.sh out=fresh_directory t=1 -Xmx512m -eoom
Checks renamed/reordered/reversed/wrapped/lowercase compressed assemblies, UTF-8
paths, duplicate-contig multiplicity, retained ambiguity symbols, empty pair lists,
circular rotations remaining unresolved, CLI parity and malformed input rejection.
Tiny synthetic files only; no biological independence or frozen split claim.
Compile complete current/ on Dori first.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaConstructionAssemblyCompareTest "$@"

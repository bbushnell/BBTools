#!/bin/bash
# Literal group closure, preserved nonmatches and complete metadata accounting.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_construction_assembly_groups_test.sh out=fresh_directory t=1 -Xmx512m -eoom
Checks inherited aliases, multiple evaluation matches, positive and unresolved
comparisons sharing an input, remaining unresolved inputs, UTF-8 provenance and
CLI parity. Malformed/omitted/duplicated/rebound metadata must reject.
Synthetic metadata only; the named FASTAs are never opened or copied.
Compile full current/ on Dori first.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaConstructionAssemblyGroupsTest "$@"

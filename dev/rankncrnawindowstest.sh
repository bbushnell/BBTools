#!/bin/bash
# Whole-CallGenes GFF/attempt parity and literal ranking-report conservation.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/rankncrnawindowstest.sh out=fresh_directory t=1 -Xmx2g -eoom
Tiny synthetic genomes and one synthetic member per native RF00174 model.
Checks full ordinary gene-path GFF/attempt parity with observation off/on,
original input namespaces, described headers, zero-window target denominators,
rank ranges and common-denominator conservation. Not an accuracy/timing panel.
Uses the current native resources and only one k/width to bound fixture work.
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
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.RankNcrnaWindowsTest "$@"

#!/bin/bash
# Compares actual short-word emissions with independent and preserved-old oracles.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	echo 'Usage: dev/hash_small_ambiguity_test.sh out=fresh.tsv expect=fixed in=optional.fa requiredefined=f t=1 -Xmx2g -eoom'
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
CP="${DIR}current/"
if [ -f "${DIR}bbtools.jar" ]; then CP="${DIR}bbtools.jar"; fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" test.cardinality.HashSmallAmbiguityTest "$@"

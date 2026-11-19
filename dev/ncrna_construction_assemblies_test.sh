#!/bin/bash
# Literal canonical-file/assembly aliases and dimension-only candidates.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_construction_assemblies_test.sh out=fresh_directory t=1 -Xmx512m -eoom
Checks explicit assembly-record versions, a real symlink alias, distinct sequences
with equal dimensions, unresolved claims and malformed original/catalog bindings.
Foreign membership paths cannot assign an assembly even with equal dimensions;
literal and canonical non-ASCII paths must round-trip as UTF-8.
Zero-locus unmeasured dimensions use the census; malformed zero combinations reject.
Uses tiny synthetic FASTAs; no training, frozen split or biological claim.
Compile complete prok/test.prok packages on Dori first.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaConstructionAssembliesTest "$@"

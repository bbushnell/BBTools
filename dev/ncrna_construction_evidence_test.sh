#!/bin/bash
# Literal construction-evidence joins and malformed provenance fixtures.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_construction_evidence_test.sh out=fresh_directory t=1 -Xmx512m -eoom
Checks cohort-qualified events, duplicate witnesses, disjoint cores, missing source
records, competing support and malformed native source/event joins.
Includes coordinate-only native support when a catalog has no saved source GFF.
Synthetic metadata only; no genome payload, training or eligibility claim.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaConstructionEvidenceTest "$@"

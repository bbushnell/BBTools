#!/bin/bash
# Literal admission, provenance conservation and assembly exclusion fixtures.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_construction_members_test.sh out=fresh_directory t=1 -Xmx512m -eoom
Tests corrected cores, duplicate witnesses, unsupported priors, seed-only fragment
retention, competing-family review, evaluation aliases, unresolved assemblies,
rare-lineage deficits, UTF-8 paths and CLI/library parity. Malformed identity,
count, diagnostic, native-event and interval inputs must fail loudly.
Synthetic metadata only; the named genome paths are never opened.
Compile the full canonical current/ tree on Dori before running.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaConstructionMembersTest "$@"

#!/bin/bash
set -eo pipefail
usage(){
	cat <<'HELP'
Usage: dev/ncrna_grade_comparison_test.sh out=new_fixture_directory

Tests fixed-reference comparison of completed grade reports using hand-counted
synthetic records. Includes input/family identity collisions, detection exchanges,
endpoint changes, masks, duplicates, exclusions and incompatible/incomplete pairs.
No caller, native search or biological reference construction runs.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.CompareNcrnaGradesTest "$@"

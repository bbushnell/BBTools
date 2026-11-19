#!/bin/bash
# Scalar counts, deterministic cross-fit, normalized weights and saved-count replay.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	echo 'Usage: dev/ncrna_partition_weight_search_test.sh out=fresh_fixture_directory t=1 -Xmx512m -eoom'
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
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaPartitionWeightSearchTest "$@"

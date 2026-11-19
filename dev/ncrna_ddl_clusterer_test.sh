#!/bin/bash
# Finite nucleotide clustering and Quantum traceback fixture, on Dori.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	echo 'Usage: dev/ncrna_ddl_clusterer_test.sh out=fresh_fixture_directory'
	exit 0
fi
test "$#" = 1
[[ $1 == out=* ]]
destination=${1#out=}
test -n "$destination"
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ];do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" test.prok.NcrnaDdlClustererTest "out=$destination"

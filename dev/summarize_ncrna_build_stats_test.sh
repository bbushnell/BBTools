#!/bin/bash
# Literal rounded-threshold uncertainty and population-conservation fixtures.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	echo 'Usage: dev/summarize_ncrna_build_stats_test.sh out=fresh_absolute_directory t=1 -Xmx256m -eoom'
	exit 0
fi
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
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" test.prok.SummarizeNcrnaBuildStatsTest "$@"

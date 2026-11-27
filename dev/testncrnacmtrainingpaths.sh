#!/usr/bin/env bash
# Independent AC paths, reverse strand, raw state anchors and disjoint ranges.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 printf '%s\n' 'Usage: dev/testncrnacmtrainingpaths.sh out=fresh_fixture_directory'
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.TraceNcrnaCmTrainingTest "$@"

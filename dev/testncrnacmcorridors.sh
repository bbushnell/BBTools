#!/usr/bin/env bash
# Literal fit/codec/fold/empty-layout fixtures, plus independent CM grammar checks.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 printf '%s\n' 'Usage: dev/testncrnacmcorridors.sh out=fresh_fixture_directory -Xmx2g -eoom'
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.BuildNcrnaCmCorridorsTest "$@"

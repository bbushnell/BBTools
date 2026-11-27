#!/usr/bin/env bash
# Literal ordinary-call census, full source/alias exclusions and model witnesses.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 printf '%s\n' 'Usage: dev/testncrnacmtraining.sh out=fresh_fixture_directory'
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.PrepareNcrnaCmTrainingTest "$@"

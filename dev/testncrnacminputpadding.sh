#!/usr/bin/env bash
# Literal resource-format, sequence-binding and sparse-fallback fixtures.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 printf '%s\n' 'Usage: dev/testncrnacminputpadding.sh out=fresh_fixture_directory t=1 -Xmx1g -eoom'
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.NcrnaCmInputPaddingTest "$@"

#!/usr/bin/env bash
# Failure-chain, cancellation and full-consumption fixtures without completion logs.
set -eo pipefail
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.prok.NcrnaChildFailureTest "installation=$DIR" "$@"

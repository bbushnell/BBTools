#!/usr/bin/env bash
# Runs the exact and exhaustive-oracle competitive selection fixtures.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CP=$DIR/current
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" prot.HbmCompetitiveGateTest

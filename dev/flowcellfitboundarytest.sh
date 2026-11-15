#!/usr/bin/env bash
# Literal degenerate-fit observations and a two-point control; no performance claim.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed -ea -eoom simd=f
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$DIR/current" \
 test.hiseq.FlowCellFitBoundaryTest "$@"

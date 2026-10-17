#!/bin/bash
# Exact guard regressions; optional in=left_right_edge.fa supplies a native proof fixture.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.GraphOverlapGateTest "$@"

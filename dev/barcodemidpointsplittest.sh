#!/bin/bash
# Literal midpoint extraction and split callers: mode=baseline|fixed case=all|api|stats|metadata.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.BarcodeMidpointSplitTest "$@"

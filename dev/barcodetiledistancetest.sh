#!/bin/bash
# Equal-core tile distance and unused legacy map row: mode=baseline|fixed case=api|legacy-extra|legacy-control.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.BarcodeTileDistanceTest "$@"

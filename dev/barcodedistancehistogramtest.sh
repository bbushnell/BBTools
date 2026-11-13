#!/bin/bash
# Finite case: mode=baseline|fixed; case=storage (default) distance=8|29|30|31 or case=report distance=3|1000.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.BarcodeDistanceHistogramTest "$@"

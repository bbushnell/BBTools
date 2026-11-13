#!/bin/bash
# Reused-string left-distance benchmark: loops=<passes> warmup=<passes> (optional).
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=128m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.BarcodeLeftDistanceBenchmark "$@"

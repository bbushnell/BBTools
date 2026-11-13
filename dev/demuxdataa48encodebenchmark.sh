#!/bin/bash
# One ordinary A48 encode timing case: length=16|27 loops=N [warmup=N].
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=128m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.DemuxDataA48EncodeBenchmark "$@"

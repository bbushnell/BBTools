#!/bin/bash
# Direct production CLI for matched barcode timing; correctness controls run separately.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=512m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" barcode.CorrelateBarcodes "$@"

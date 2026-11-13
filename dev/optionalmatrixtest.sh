#!/bin/bash
# Optional barcode implementation checks; no private implementation is shipped here.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/${BARCODE_OPTIONAL_CLASSPATH:+:$BARCODE_OPTIONAL_CLASSPATH}" test.barcode.OptionalMatrixTest "$@"

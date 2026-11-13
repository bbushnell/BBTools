#!/bin/bash
# One finite CountBarcodes2 client-fallback case: mode=baseline|fixed case=<name>.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.CountBarcodes2UseServerTest "$@"

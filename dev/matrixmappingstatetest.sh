#!/bin/bash
# One finite mapping-state case: mode=baseline|fixed driver=serial|local|shared seeded=t|f.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.MatrixMappingStateTest "$@"

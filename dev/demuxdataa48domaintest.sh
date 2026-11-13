#!/bin/bash
# Finite base-5 domain checks: mode=baseline|fixed. Assertions are required.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.DemuxDataA48DomainTest "$@"

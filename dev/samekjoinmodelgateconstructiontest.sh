#!/bin/bash
# One fresh Linux JVM per model-gate constructor case; outdir must be new.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.SameKJoinModelGateConstructionTest "$@"

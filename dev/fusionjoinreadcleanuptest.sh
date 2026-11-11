#!/bin/bash
# One fresh Linux JVM per trace-reader cleanup case; outdir must not exist.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=512m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.FusionJoinReadCleanupTest "$@"

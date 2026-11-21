#!/usr/bin/env bash
# Actual processor header conversion; run fresh matched JVMs with a measured iteration count.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=256m --mode=fixed simd=f "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.processor.UnderscoreToSpaceBenchmark "$@"

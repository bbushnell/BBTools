#!/usr/bin/env bash
# Actual CountGC reader parsing; use matched fresh JVMs and a calibrated pass count.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=256m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.CountGCParsingBenchmark "$@"

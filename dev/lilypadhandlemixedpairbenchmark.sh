#!/bin/bash
# Matched Lilypad mixed-pair edge update benchmark.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$DIR/current/" \
	test.consensus.LilypadHandleMixedPairBenchmark "$@"

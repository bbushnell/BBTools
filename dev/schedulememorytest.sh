#!/bin/bash
# Scalar schedule regression; no kmer tables or read corpus are allocated.
# javasetup probes optional environment variables, as in other BBTools wrappers.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
CP="${SCHEDULE_MEMORY_TEST_CLASSPATH:-$DIR/current/}"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" kmer.ScheduleMakerMemoryTest "$@"

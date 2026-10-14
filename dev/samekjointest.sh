#!/bin/bash
# Deterministic reference-oracle and candidate-geometry fixtures.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
CP="${SAMEK_TEST_CLASSPATH:-$DIR/current/}"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" assemble.SameKJoinTest "$@"

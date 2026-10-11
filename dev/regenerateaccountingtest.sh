#!/bin/bash
# Small deterministic explicit-table accounting tests; no read corpus required.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
CP="${REGENERATE_TEST_CLASSPATH:-$DIR/current/}"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" ukmer.HashArrayU1DRegenerateTest "$@"

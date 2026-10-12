#!/bin/bash
# Unit tests use fake child launchers; no real assembly workload runs locally.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
CP="${TADPIPE_TEST_CLASSPATH:-$DIR/current/}"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" assemble.TadPipeTest "$@"

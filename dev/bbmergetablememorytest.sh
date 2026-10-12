#!/bin/bash
# Constructor-only tests; no table loading or assembly work occurs locally.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
CP="${BBMERGE_MEMORY_TEST_CLASSPATH:-$DIR/current/}"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" jgi.BBMergeTableMemoryTest "$@"

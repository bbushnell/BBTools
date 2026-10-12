#!/bin/bash
# Scalar sketch tests and deterministic fixture generation; no heavy local run.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
CP="${KHIST_TEST_CLASSPATH:-$DIR/current}"
. "$DIR/javasetup.sh"
. "$DIR/memdetect.sh"
parseJavaArgs "--xmx=256m" "--xms=64m" "--mode=fixed" "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$CP" cardinality.KmerHistMultiTest "$@"

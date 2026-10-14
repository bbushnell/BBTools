#!/bin/bash
# Tiny schema fixtures, not simulated-genome truth or model-performance tests.
# generate new_root [variant] ; check converted_output
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" assemble.SameKTrainingCorpusTest "$@"

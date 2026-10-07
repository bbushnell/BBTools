#!/bin/bash
# Run after a full dna/gff/prok compile; fixture directory is required and caller-owned.
set -eo pipefail
DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CP="${BBTOOLS_TEST_CLASSES:-$DIR/current}"
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" prok.CallGenesGeneticCodeTest "$@"

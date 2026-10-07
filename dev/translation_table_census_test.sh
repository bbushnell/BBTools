#!/bin/bash
# Uses freshly compiled production classes and writes fixtures outside the code tree.
set -eo pipefail
DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CP="${BBTOOLS_TEST_CLASSES:-$DIR/current}"
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" prok.TranslationTableCensusTest "$@"

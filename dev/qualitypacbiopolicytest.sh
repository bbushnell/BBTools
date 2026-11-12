#!/bin/bash
# One quality-histogram, PacBio helper-policy, or PacBio helper benchmark check.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=64m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.align2.QualityPacBioPolicyTest "$@"

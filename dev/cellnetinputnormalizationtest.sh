#!/bin/bash
# Analytic parser/copy/round-trip checks for native input standardization.
set -o pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}") || exit 1
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd) || exit 1
if [[ -f "$DIR/bbtools.jar" ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current/; fi
source "$DIR/javasetup.sh" || exit 1
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@" || exit 1
setEnvironment || exit 1
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" ml.CellNetInputNormalizationTest

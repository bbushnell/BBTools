#!/bin/bash
# One finite output-policy case per fresh JVM; dir must not already exist.
# mode=baseline|fixed report=summary|labels scenario=fresh|reject|replace|append dir=NEW
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.CompareLabelsOutputPolicyTest "$@"

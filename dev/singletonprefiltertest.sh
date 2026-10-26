#!/bin/bash
# Native Solid64 singleton filtering: exact predicate parity under crowded filters.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.ukmer.SingletonPrefilterTest "$1"

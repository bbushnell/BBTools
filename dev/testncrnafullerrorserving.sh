#!/usr/bin/env bash
# Fixed-vector package/refiner parity using existing exports; never fits a net.
set -eo pipefail
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.NcrnaFullErrorServingTest "$@"

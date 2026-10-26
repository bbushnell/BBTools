#!/bin/bash
# Endpoint lookup parity against packed Hollow and alternate-span diagnostic queries.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.ukmer.HollowKmerLookupTest

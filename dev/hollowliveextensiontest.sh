#!/bin/bash
# Actual Tadpole2 entry-point parity with controlled nominal and alternate counts.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.HollowLiveExtensionTest in="$DIR/testdata/hollow_rollback_queries.fa"

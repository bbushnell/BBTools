#!/bin/bash
# Controlled in-memory stop/recovery/rollback transactions; no new read files.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.HollowAnchorRescueTest in="$DIR/testdata/hollow_rollback_queries.fa"

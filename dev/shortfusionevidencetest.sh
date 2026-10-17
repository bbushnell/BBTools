#!/bin/bash
# Rejection-only short-end evidence regressions; assertions always enabled.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "${SHORT_FUSION_TEST_CLASSPATH:-$DIR/current/}" assemble.ShortFusionEvidenceTest "$@"

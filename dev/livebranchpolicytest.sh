#!/bin/bash
# Deterministic policy/feature invariance tests; no assembly behavior is changed.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
CP="${LIVE_BRANCH_TEST_CLASSPATH:-$DIR/current/}"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" assemble.LiveBranchPolicyTest "$@"

#!/bin/bash
# Synthetic final-only branch veto tests; optional assembly args test policy restoration.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "${GRAPH_BRANCH_TEST_CLASSPATH:-$DIR/current/}" assemble.GraphBranchGateTest "$@"

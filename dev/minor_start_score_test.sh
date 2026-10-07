#!/usr/bin/env bash
# Runs after the full affected-package build. Set JAVA_TOOL_OPTIONS to select the JVM policy.
set -eo pipefail
DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CP=${BBTOOLS_TEST_CLASSES:-$DIR/current}
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" prok.MinorStartScoreTest "$@"

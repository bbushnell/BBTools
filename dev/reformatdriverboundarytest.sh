#!/usr/bin/env bash
# Finite Reformat driver contracts; each invocation uses a fresh JVM and directory.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.ReformatDriverBoundaryTest "$@"

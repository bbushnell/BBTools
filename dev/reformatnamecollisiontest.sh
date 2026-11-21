#!/usr/bin/env bash
# Literal name allocation/scope and real-driver trimming controls; use a new outdir.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.processor.ReformatNameCollisionTest "$@"

#!/usr/bin/env bash
# Ordinary two-worker BBMask caller controls; no production compilation.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed -ea simd=f "$@"
setEnvironment
exec java -Dfile.encoding=UTF-8 -XX:ActiveProcessorCount=2 $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.BBMaskWorkerJoinControlsTest "bbtools=$DIR" "$@"

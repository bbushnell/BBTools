#!/usr/bin/env bash
# Owned filesystem-failure baselines and a literal concatenation control; no compilation.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed -ea simd=f "$@"
setEnvironment
exec java -Dfile.encoding=UTF-8 -XX:ActiveProcessorCount=1 $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.BBMergeMultiTempCleanupTest "$@"

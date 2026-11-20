#!/usr/bin/env bash
# Two AdjustHomopolymers processing exceptional-cleanup children; no compilation or deployment.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
exec java -Dfile.encoding=UTF-8 -XX:ActiveProcessorCount=1 $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.AdjustHomopolymersProcessingCleanupTest "mode=parent" "tree=$DIR" "$@"

#!/bin/bash
# Exact conditioned-distribution oracle; no genomes or network training.
set -eo pipefail
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$DIR/current/" test.prok.NcrnaBoundaryJitterTest "$@"

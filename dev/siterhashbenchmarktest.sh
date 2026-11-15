#!/usr/bin/env bash
# SiteR key hash controls and caller-shaped HashMap lookup benchmark.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed -ea -eoom simd=f
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 \
	-cp "$DIR/current" test.pacbio.SiteRHashBenchmarkTest "$@"

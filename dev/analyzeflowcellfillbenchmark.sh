#!/usr/bin/env bash
# CPU-stage workload; calibrate to >=10 seconds before drawing timing conclusions.
set -eo pipefail
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=256m --mode=fixed -ea -eoom simd=f
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS -XX:ActiveProcessorCount=1 -cp "$DIR/current" \
 test.hiseq.AnalyzeFlowCellFillBenchmark "$@"

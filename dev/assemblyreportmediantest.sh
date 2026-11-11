#!/bin/bash
# Seven direct report-median cases; output uses exact hexadecimal double bits.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=64m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.AssemblyReportMedianTest "$@"

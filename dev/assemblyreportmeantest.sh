#!/bin/bash
# Direct timing-mean contracts: nine valid raw-bit cases and nine invalid inputs.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=64m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.assemble.AssemblyReportMeanTest "$@"

#!/usr/bin/env bash
# Run the finite CalcTrueQuality writeMatrix close-status fixture.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=64m --xms=32m --mode=fixed simd=f eoom "$@"
setEnvironment
exec java -Dfile.encoding=UTF-8 -XX:ActiveProcessorCount=1 $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.CalcTrueQualityWriteMatrixCloseTest "$@"

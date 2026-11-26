#!/usr/bin/env bash
# Generate a scaled CalcTrueQuality q12 timing SAM plus independent matrix oracles.
set -eo pipefail

DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"

parseJavaArgs --xmx=256m --xms=32m --mode=fixed simd=f eoom "$@"
setEnvironment
exec java -Dfile.encoding=UTF-8 -XX:ActiveProcessorCount=1 $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.CalcTrueQualityQ12TimingData "$@"

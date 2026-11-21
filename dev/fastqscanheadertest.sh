#!/usr/bin/env bash
# Literal raw FASTQ header/trim controls; run from the BBTools root with a new outdir.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -Dfile.encoding=UTF-8 -cp "$DIR/current/" test.stream.FastqScanHeaderTest "$@"

#!/bin/bash
# Deterministic NovaDemux subset timing FASTQ and legacy-report verification.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.barcode.NovaSubsetBenchmarkInput "$@"

#!/usr/bin/env bash
# Independent fixed identifier-scope benchmark specification; writes no FASTQ.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.processor.ReformatNamingOutputOracle "$@"

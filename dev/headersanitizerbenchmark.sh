#!/usr/bin/env bash
# Actual legacy header sanitizers; short invocations are calibration, not timing evidence.
set -eo pipefail
if [[ $# = 0 || $1 = -h || $1 = --help ]]; then
 cat <<'USAGE'
Usage: headersanitizerbenchmark.sh phase=baseline|fixed target=legacy|legacyold profile=allowed|dirty iterations=1m
Use JDK_JAVA_OPTIONS='-XX:ActiveProcessorCount=1' for each fresh single-CPU JVM.
Each profile has64 distinct128-char ASCII headers; matched final samples need >=10CPU seconds.
USAGE
 exit 0
fi
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=256m --mode=fixed simd=f "$@"
setEnvironment
java -Dfile.encoding=UTF-8 $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.processor.HeaderSanitizerBenchmark "$@"

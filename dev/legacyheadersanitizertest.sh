#!/usr/bin/env bash
# Direct modern and legacy String sanitizer controls with explicit UTF-8.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java -Dfile.encoding=UTF-8 $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.LegacyHeaderSanitizerTest "$@"

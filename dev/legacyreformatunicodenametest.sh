#!/usr/bin/env bash
# Direct legacy String Unicode controls with an explicit UTF-8 baseline.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=32m --mode=fixed "$@"
setEnvironment
java -Dfile.encoding=UTF-8 $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.LegacyReformatUnicodeNameTest "$@"

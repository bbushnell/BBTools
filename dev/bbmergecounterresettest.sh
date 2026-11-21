#!/usr/bin/env bash
# Two real BBMerge phases in one fresh JVM, followed by one direct reset control.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=128m --xms=32m --mode=fixed -ea simd=f "$@"
setEnvironment
exec java -Dfile.encoding=UTF-8 -Duser.language=en -Duser.country=US -XX:ActiveProcessorCount=1 $EA $EOOM $SIMD $XMX $XMS -cp "$DIR/current/" test.jgi.BBMergeCounterResetTest "$@"

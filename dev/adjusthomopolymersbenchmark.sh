#!/usr/bin/env bash
# One fresh-JVM pilot/sample against retained compiled production classes; no compilation.
set -eo pipefail
DIR="$(cd "$(dirname "$0")/.." && pwd)"
TARGET=""
for arg in "$@"; do
 if [[ "$arg" == targetclasses=* ]]; then TARGET="${arg#targetclasses=}"; fi
done
[[ -n "$TARGET" && -d "$TARGET" && -f "$TARGET/jgi/AdjustHomopolymers.class" ]]
. "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=256m --mode=fixed simd=f "$@"
setEnvironment
exec java -Dfile.encoding=UTF-8 -XX:ActiveProcessorCount=1 $EA $EOOM $SIMD $XMX $XMS -cp "$TARGET:$DIR/current/" test.jgi.AdjustHomopolymersBenchmark "$@"

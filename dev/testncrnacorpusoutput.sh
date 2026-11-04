#!/usr/bin/env bash
set -eo pipefail
class=test.prok.NcrnaCorpusOutputTest
if [[ "${1:-}" = log ]];then class=test.prok.NcrnaCorpusLogTest;shift;fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed -eoom "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$class" "$@"

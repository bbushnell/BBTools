#!/usr/bin/env bash
set -eo pipefail
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
main=test.prok.NcrnaSurveyOutputTest
if [[ "${1:-}" = workers ]];then main=test.prok.NcrnaSurveyWorkerTest;shift;fi
if [[ "${1:-}" = cli ]];then main=test.prok.NcrnaSurveyCliTest;shift;fi
parseJavaArgs --xmx=512m --xms=64m --mode=fixed -eoom "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$main" "$@"

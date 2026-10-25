#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=2g --mode=fixed "$@";then exit 1;fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
exec java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -cp "$root/current" test.prok.NcrnaMitoConsensusTest "$@"

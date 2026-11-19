#!/usr/bin/env bash
# Diagnostic comparison only; no network evaluation or acceptance tolerance.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 printf '%s\n' 'Usage: comparencrnamlcaptures.sh in=before.tsv in2=after.tsv out=fresh_report.tsv ow=f'
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
ARGS=()
for arg in "$@";do
 case "$arg" in
  -Xmx*|-Xms*|-ea|-da|ea|da|-eoom|eoom) ;;
  *) ARGS+=("$arg") ;;
 esac
done
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$DIR/current" test.prok.NcrnaMlCaptureCompare "${ARGS[@]}"

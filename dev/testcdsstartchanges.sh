#!/usr/bin/env bash
# Five-input finite source-relative CDS grading fixture; no external genomes.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 printf '%s\n' 'Usage: dev/testcdsstartchanges.sh out=fresh_directory t=2 -Xmx512m'
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=512m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.CdsStartChangesTest "$@"

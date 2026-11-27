#!/usr/bin/env bash
# Aggregate exact training paths and glocal anchors; never calls genes or fits.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 printf '%s\n' 'Usage: dev/tracencrnacmtraining.sh in=training.calls.tsv.gz bindings=training.sources.tsv refs=main.fa,second.fa,third.fa maps=main.tsv,second.tsv models=main.cm,second.cm expectedcalls=N expectedinputs=N first=0 limit=64 padding=20 maxcells=500m plan=f out=fresh_directory t=1 -Xmx8g'
 printf '%s\n' 'The pilot owns its explicit call-index range; later jobs use disjoint ranges. plan=t reports matrix bounds without searching. Pass trd=f itn=f utot=f dotdashxton=f lctn=f to preserve native records.'
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=8g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.TraceNcrnaCmTraining "$@"

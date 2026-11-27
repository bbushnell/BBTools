#!/usr/bin/env bash
# Census and source-qualified inputs from completed ordinary calls; no alignment.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 printf '%s\n' 'Usage: dev/preparencrnacmtraining.sh in=prefixes.tsv bindings=full_sources.tsv audits=main,second,third exclude=evaluation.tsv refs=main.fa,second.fa,third.fa expectedmodels=412 expectedinputs=40686 expectedprefixes=32 mincalls=100 minassemblies=10 out=fresh_prefix'
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.PrepareNcrnaCmTraining "$@"

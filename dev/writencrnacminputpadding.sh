#!/usr/bin/env bash
# Export fixed padding pairs from one consensus FASTA, without fitting support.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 printf '%s\n' 'Usage: dev/writencrnacminputpadding.sh in=consensus.fa out=fresh_directory pad5=20,25,30 pad3=20,10,10 t=1 -Xmx1g -eoom'
 printf '%s\n' 'Writes padding_FIVE_THREE.tsv for each parallel pair. Every current model uses the constants, including models without teaching support. Exact model order, names and consensus bases bind each resource. No teaching data, fitting, caller or search.'
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.WriteNcrnaCmInputPadding "$@"

#!/usr/bin/env bash
# Reports complete three-arm family transitions, not biological accuracy.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 echo 'Usage: hbmexpansiontransitions.sh config=FILE out=NEW_DIRECTORY expected=N oldfamilies=N'
 echo 'Required inputs, each paired with KEYsha80=PIN: identity, shipping, restricted, expanded, map.'
 echo 'Inputs have matching comparison_gN IDs in original order; map names bacteria_0..49/archaea_0..49.'
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
source "$DIR/javasetup.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
 case ${arg,,} in
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$DIR/current" prot.HbmExpansionTransitions "${ARGS[@]}"

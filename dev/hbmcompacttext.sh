#!/usr/bin/env bash
# Read compact HBM text and verify native graph and optional-header parity.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
  echo 'Usage: hbmcompacttext.sh in=model.hbmc(.gz) [baseline=accepted.hbmt.gz provenance=manifest.tsv] [profile=schema7.tsv.gz]'
  echo 'selftest=t tests metadata omission, topology, counts and malformed input.'
  exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -f $DIR/bbtools.jar ]]; then CP=$DIR/bbtools.jar; else CP=$DIR/current; fi
source "$DIR/javasetup.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
  case ${arg,,} in
    --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
    -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
    simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
    *=*) key=${arg%%=*}; ARGS+=("${key,,}=${arg#*=}");;
    *) ARGS+=("$arg");;
  esac
done
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" prot.HbmCompactTextProbe "${ARGS[@]}"

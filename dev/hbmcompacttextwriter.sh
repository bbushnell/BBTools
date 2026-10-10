#!/usr/bin/env bash
# Write experimental compact HBM text without changing runtime models.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
  echo 'Usage: hbmcompacttextwriter.sh in=accepted.hbmt.gz out=NEW.txt order=frequency|rows|alphabetical'
  echo 'Optional filter=none|singletons|states|both stats=NEW.tsv; selftest=t runs filter boundaries.'
  echo 'Optional profile=schema7.tsv(.gz) imports effective cutoffs and length bounds; output should use .hbmc(.gz).'
  echo 'Writer-only size experiment; alphabets are computed after filtering. No reader or accuracy claim.'
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
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" prot.HbmCompactTextWriter "${ARGS[@]}"

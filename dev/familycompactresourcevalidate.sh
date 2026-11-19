#!/usr/bin/env bash
# Validate a rendered compact-family resource against its source references and sidecars.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
  echo 'Usage: familycompactresourcevalidate.sh mode=compare view=short references=old.tsv.gz compact=families.tsv.gz profilevectors=profiles.tsv.gz lengthmedians=family_lengths.tsv familyhydro=member_hydrophobicity.tsv familygc=family_gc.tsv familygcmap=family_map.tsv scalarbounds=scalar_bounds.tsv alphabet=alphabet.config matrixbinding=POSSUM55:... out=NEW_DIR [t=64]'
  echo '       familycompactresourcevalidate.sh mode=blacklist short=shared.faa tails=full.faa out=blacklist.faa'
  echo '       familycompactresourcevalidate.sh mode=matrixbinding matrixresource=resources/bbgap/POSSUM55.tsv [acceptedmatrix=25to55.tsv]'
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
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" prot.FamilyCompactResourceValidate "${ARGS[@]}"

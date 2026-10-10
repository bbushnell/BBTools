#!/usr/bin/env bash
# Builds an experimental profile from one pinned, complete seed family.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 cat <<'HELP'
Usage: newfamilyhbm.sh manifest=family_manifest.tsv manifestsha80=PIN source=DIR
  rank=N background=background.tsv backgroundsha80=PIN runtime=RUNTIME out=NEW_DIR
All family members are used. BLOSUM62 seed placement precedes two positional
profile passes (beta .01, clipping -4..11, gap4). Padding starts at20 and grows
when required. Inputs and native bundle round trips are verified. No installation.
HELP
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -d $DIR/current ]]; then CP=$DIR/current; else CP=$DIR/bbtools.jar; fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
 case ${arg,,} in
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
parseJavaArgs --xmx=4g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" prot.NewFamilyHbmBuilder "${ARGS[@]}"

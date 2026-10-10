#!/usr/bin/env bash
# Refines changed nonempty families and preserves unchanged or unresolved-empty seeds.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 echo 'Usage: hbmreassignedprofiles.sh old=MANIFEST oldsha80=PIN new=MANIFEST newsha80=PIN models=REGISTRY modelssha80=PIN background=FILE backgroundsha80=PIN runtime=DIR [ranks=FILE] out=NEW_DIRECTORY'
 echo 'Uses the frozen original background; EMPTY_REQUIRES_DECISION prevents silent model deletion.'
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CP=$DIR/current
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
parseJavaArgs --xmx=8g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" prot.HbmReassignedProfiles "${ARGS[@]}"

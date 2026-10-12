#!/usr/bin/env bash
# Verifies and combines old refinements with new experimental seed profiles.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 cat <<'HELP'
Usage: hbmfamilyunion.sh resources=OLD_PRODUCTION_RESOURCES
  oldmanifest=FILE oldmanifestsha80=PIN oldsource=DIR oldruntime=DIR oldfamilies=DIR
  newmanifest=FILE newmanifestsha80=PIN newsource=DIR newruntime=DIR newfamilies=DIR
  background=FILE backgroundsha80=PIN runtime=THIS_RUNTIME out=NEW_DIRECTORY
Optional oldranks=FILE and newranks=FILE select zero-based row offsets within
their respective manifests for bounded checks. Defaults include every family.
The output keeps all selected graphs unchanged, retains permanent IDs, emits
dense consensus indexes for the sidecar, and verifies a native bundle round trip.
The position-score background remains the explicitly pinned original background.
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
parseJavaArgs --xmx=24g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" prot.HbmFamilyUnion "${ARGS[@]}"

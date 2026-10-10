#!/usr/bin/env bash
# Recounts all completed compressed assignments against their bound model metadata.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 echo 'Usage: hbmassignmentcollect.sh queries=FILE queriessha80=PIN modelconfig=FILE modelconfigsha80=PIN root=FLEET_ROOT expected=N out=NEW_DIRECTORY'
 echo 'Writes query_ids.txt for a separate global uniqueness check; PASS explicitly leaves that check pending.'
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CP=$DIR/current
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=(); main=prot.HbmAssignmentCollect
for arg in "$@"; do
 case ${arg,,} in
  selftest=t) main=prot.HbmAssignmentCollectTest;;
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
parseJavaArgs --xmx=2g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" "$main" "${ARGS[@]}"

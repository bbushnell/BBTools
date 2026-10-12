#!/usr/bin/env bash
# Compares family memberships exactly by ID and encoded residues.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 echo 'Usage: hbmmembershipcompare.sh old=MANIFEST oldsha80=PIN new=MANIFEST newsha80=PIN [ranks=FILE] out=NEW_DIRECTORY'
 echo 'Reports UNCHANGED, CHANGED or EMPTY; zero-member models are not automatically discarded.'
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
CP=$DIR/current
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=(); main=prot.HbmMembershipCompare
for arg in "$@"; do
 case ${arg,,} in
  selftest=t) main=prot.HbmMembershipCompareTest;;
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
parseJavaArgs --xmx=2g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" "$main" "${ARGS[@]}"

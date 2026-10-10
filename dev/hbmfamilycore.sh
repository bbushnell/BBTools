#!/usr/bin/env bash
# Recomputes paired half-member cores with the final positional profiles.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]]; then
 cat <<'HELP'
Usage: hbmfamilycore.sh manifest=FILE manifestsha80=PIN models=FILE modelssha80=PIN
  consensus=FILE consensussha80=PIN background=FILE backgroundsha80=PIN
  ranks=FILE out=NEW_DIRECTORY
The member manifest contains absolute FASTA paths. Models is the verified union's
family_artifacts.tsv. Ranks contains dense source indexes, one per line.
Every member is aligned with logodds profiles, beta .01, clipping -4..11, gap4,
and the explicitly pinned background. Only paired columns increment depth.
The core spans the first through last columns covered by at least half the members.
Use selftest=t for the exact paired-path and endpoint fixtures.
HELP
 exit 0
fi
SCRIPT=$(readlink -f -- "${BASH_SOURCE[0]}")
DIR=$(cd "$(dirname "$SCRIPT")/.." && pwd)
if [[ -d $DIR/current ]]; then CP=$DIR/current; else CP=$DIR/bbtools.jar; fi
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=(); main=prot.HbmFamilyCore
for arg in "$@"; do
 case ${arg,,} in
  selftest=t) main=prot.HbmFamilyCoreTest;;
  --xmx=*|--xms=*|--mem=*|--percent=*|--mode=*|-xmx*|xmx=*|-xms*|xms=*|\
  -ea|-da|ea|da|-eoom|eoom|exitonoutofmemoryerror|-exitonoutofmemoryerror|\
  simd|simd=t|simd=true|simd=f|simd=false|nosimd) JVM_ARGS+=("$arg");;
  *) ARGS+=("$arg");;
 esac
done
parseJavaArgs --xmx=4g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -XX:ActiveProcessorCount=1 -cp "$CP" "$main" "${ARGS[@]}"

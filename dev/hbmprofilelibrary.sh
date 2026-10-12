#!/usr/bin/env bash
# Offline sharded profile refinement; no production installation.
set -eo pipefail
if [[ $# == 0 ]]; then
 echo 'Usage: hbmprofilelibrary.sh mode=plan|shard|combine resources=DIR manifest=FILE manifestsha80=PIN source=FASTA_DIR out=NEW_DIR'
 echo 'plan: shards=32; shard: ranks=FILE; combine: families=DIR [ranks=FILE for explicit subset]'
 exit 0
fi
DIR=$(cd "$(dirname "$0")/.." && pwd)
CP=$DIR/current
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=(); main=prot.HbmProfileLibrary
for arg in "$@"; do
 case "$arg" in -Xmx*|-Xms*|-ea|-da|-eoom) JVM_ARGS+=("$arg");; selftest=t) main=prot.HbmProfileLibraryTest;; *) ARGS+=("$arg");; esac
done
parseJavaArgs --xmx=24g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" "$main" "${ARGS[@]}" "runtime=$DIR"

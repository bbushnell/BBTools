#!/usr/bin/env bash
# Experimental single-family profile refinement; writes fresh outputs only.
set -eo pipefail
if [[ $# == 0 ]]; then
 echo 'Usage: hbmprofilepilot.sh in=family.faa resources=DIR manifest=FILE manifestsha80=PIN rank=N out=NEW_DIR padding=20'
 exit 0
fi
DIR=$(cd "$(dirname "$0")/.." && pwd)
CP=$DIR/current
source "$DIR/javasetup.sh"
source "$DIR/memdetect.sh"
ARGS=(); JVM_ARGS=()
for arg in "$@"; do
 case "$arg" in -Xmx*|-Xms*|-ea|-da|-eoom) JVM_ARGS+=("$arg");; *) ARGS+=("$arg");; esac
done
parseJavaArgs --xmx=16g --xms=256m --mode=fixed "${JVM_ARGS[@]}"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $PROXY -cp "$CP" prot.HbmProfilePilot "${ARGS[@]}" "runtime=$DIR"

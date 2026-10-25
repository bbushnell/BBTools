#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Exhaustive native best-model identity; no seeding, shortlist or early acceptance.'
  echo 'in=MEMBERS.fa.gz ref=MODELS.fa.gz tree=TREE outprefix=NEW_PREFIX arm=NAME expected=COPIES expectedmodels=N verified=INPUTS.sha80 t=16 -Xmx8g'
  exit 0
fi
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=8g --mode=fixed "$@";then exit 1;fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
exec java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=16 -cp "$root/current" prok.NcrnaMitoCrossLibrary "$@"

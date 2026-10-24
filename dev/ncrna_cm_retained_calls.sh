#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Compare complete Infernal1.1.5 tblout raw multisets, ignoring only E-value changes.'
  echo 'reference=TABLE candidate=TABLE model=RFNNNNN clen=N ga=N outprefix=NEW_PREFIX -Xmx1g'
  echo 'All other fixed fields and duplicate multiplicity must match. Changed rows are retained, not silently paired by overlap.'
  exit 0
fi
DIR=$root
set +u
source "$DIR/javasetup.sh"
if ! parseJavaArgs --xmx=1g --mode=fixed "$@";then exit 1;fi
setEnvironment
set -u
read -r -a java_options <<< "$EA $EOOM $SIMD $XMX $XMS $XSS"
exec java "${java_options[@]}" -XX:+ExitOnOutOfMemoryError -XX:ActiveProcessorCount=2 -cp "$root/current" prok.NcrnaCmRetainedCalls "$@"

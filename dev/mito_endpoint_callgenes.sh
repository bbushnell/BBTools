#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if [[ $# = 0 || "${1:-}" = --help ]];then
  echo 'Experimental actual CallGenes: normal flags plus mitonet=FILE mitotables=FILE mitolibrarysha80=PIN mitoradius=192.'
  echo 'Omit all mito resource arguments for NN-off control. Only added mitochondrial models are refined; original five16S models and other families retain their paths.'
  exit 0
fi
source "$root/calcmem.sh"
z=-Xmx8g
# The standard environment helper deliberately probes optional unset variables.
set +eu
setEnvironment
parseXmx "$@"
set -eu
vectorFlags=()
if [[ -n "$SIMD" ]];then read -r -a vectorFlags <<< "$SIMD";fi
exec java -ea -XX:+ExitOnOutOfMemoryError "$z" "${vectorFlags[@]}" -cp "$root/current" test.prok.MitoEndpointCallGenes "$@"

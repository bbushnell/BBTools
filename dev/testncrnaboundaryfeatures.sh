#!/usr/bin/env bash
# Runs the existing literal feature-layout controls with native JVM setup.
set -eo pipefail
if [[ $# -eq 0 || "$1" == -h || "$1" == --help ]]; then
  echo 'Usage: bash dev/testncrnaboundaryfeatures.sh v2|v3 [JVM options]'
  echo 'Checks feature layout and absence of per-candidate alignment on finite inputs.'
  exit 0
fi
case "$1" in
  v2) main=test.prok.NcrnaBoundaryFeatureV2Test ;;
  v3) main=test.prok.NcrnaBoundaryFeatureV3Test ;;
  *) echo 'Expected v2 or v3.' >&2; exit 1 ;;
esac
shift
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$main" "$@"

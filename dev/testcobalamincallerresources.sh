#!/usr/bin/env bash
# Finite actual-data package/lazy-worker checks; no genome benchmark.
set -eo pipefail
if [[ $# == 0 || ${1:-} == --help ]];then
 printf '%s\n' 'Usage: dev/testcobalamincallerresources.sh in=serving_group out=fresh_fixture_directory scorecell=padding|greedy' 'Compares ordinary lazy setup with explicit development setup on a discovery-positive consensus on both strands. The actual package must contain zero (padding) or ten (greedy) relative model-factor changes; the greedy probe uses a changed main-member coefficient. Checks disabled loading, missing/mismatched group roles, main covariance selection before old tips, and explicit endpoint conflicts. No genome accuracy or timing claim.'
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.CobalaminCallerResourcesTest "$@"

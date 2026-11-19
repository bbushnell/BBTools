#!/usr/bin/env bash
# Saved-table clan pooling only; no calling, inference, fitting or selection.
set -eo pipefail
class=prok.NcrnaCenteredGatePool
if [[ ${1:-} == test ]];then
 class=test.prok.NcrnaCenteredGatePoolTest
 shift
fi
if [[ $# == 0 || ${1:-} == -h || ${1:-} == --help ]];then
 printf '%s\n' 'Usage: ncrnacenteredgatepool.sh in=summaries.tsv in2=witnesses.tsv out=fresh_directory' \
  'Or: ncrnacenteredgatepool.sh test out=fresh_fixture_directory' \
  'Pools exact saved counts and square sums with unchanged RF01689; retains every tied cutoff combination.'
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" "$class" "$@"

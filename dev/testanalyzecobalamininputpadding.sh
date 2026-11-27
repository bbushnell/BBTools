#!/usr/bin/env bash
# Tiny literal saved-evidence fixtures; no caller, CM search or real genome input.
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/testanalyzecobalamininputpadding.sh out=fresh_fixture_directory
Tests fixed reference intersections, alias-group separation, both strands,
model switches, exact squared-error ties, clipping and search-work arithmetic.
Invalid alias identity, witness offsets, missing contig terminal, profile counts,
consensus sequence and duplicate reference assignments must fail.
No real caller, covariance search, grading or scientific measurement runs.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.AnalyzeCobalaminInputPaddingTest "$@"

#!/usr/bin/env bash
# Finite saved-proposal/reference oracles; no caller, CM search or genome input.
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/testanalyzecobalamincmproposals.sh out=fresh_fixture_directory
Tests proposal analysis, the explicit-member score-guard report, the reference
codec and saved acceptance curves with separate source-group splits.
Literal fixtures include both strands, exact endpoints protected by the guard,
inexact endpoints a rejected proposal would fix, selected and unselected calls,
masked/ambiguous reference oracles, untraced members and missing profile counts.
Acceptance fixtures distinguish verifier and endpoint identities, inclusive ties,
matched-inexact subsets, measured-CM identity versus protected assignments, native
score gaps, unrelated reference families and missing-score accounting.
No caller, covariance search, real genome or scientific measurement is run.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" test.prok.AnalyzeCobalaminCmProposalsTest "$@"

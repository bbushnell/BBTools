#!/usr/bin/env bash
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/summarizencrnaqueryoutcomes.sh in=completed_grade_prefix familymap=original_map.tsv family=cobalamin expected=1859,425,39 out=fresh_summary.tsv
Counts recorded joint GradeNcrnaReference query outcomes by emitted model prefix.
Included calls have four disjoint buckets: TP, MASKED, unmatched duplicate and
other unmatched. Independent duplicate/scored-eligibility flags remain counters.
Only included prefixes of the selected scoring family are reported, including
zero-query prefixes. Excluded/reference rows never enter their call totals.
expected= optionally checks included counts in original map order.
All counts must reconcile with the original
joint family summary before output opens. No matching, sequence/GFF read or
biological subtype assignment occurs. Text retains the grade's percent escaping.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.SummarizeNcrnaQueryOutcomes "$@"

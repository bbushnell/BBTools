#!/usr/bin/env bash
# Default: literal join fixtures. Optional report arguments run the diagnostic only.
# bash dev/bbgap_callgenes_fate_report_test.sh BUILD mode=changes manifest=INPUT out=PREFIX
# Changes manifest: genome<TAB>fate<TAB>comparison; fate is the candidate arm's
# support-fate-v1 TSV, comparison is ProteinStopEvaluation's stops-v1 .stops.tsv.
# Outputs: PREFIX.changes.tsv and PREFIX.summary.tsv. Keys, not reference rows,
# are counted. Untraced means no supported-stop trace, not an inferred loss cause.
set -euo pipefail
build=${1:?frozen build}
shift
[[ $build == /tmp/bbgap-callgenes-* && -f $build/PASS ]]
target=test.prot.ProteinSupportFateReportTest
if (( $# )); then target=prot.ProteinSupportFateReport; fi
exec nice -n10 java -ea -Xmx512m -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector \
 -cp "$build/runtime/current" "$target" "$@"

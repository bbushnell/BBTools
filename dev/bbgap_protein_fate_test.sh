#!/usr/bin/env bash
# Tests observational fate accounting and early CLI rejection without loading models.
set -euo pipefail
build=${1:?immutable full CallGenes build}
[[ $build == /tmp/bbgap-callgenes-* && -f $build/PASS ]]
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
bash "$repo/dev/bbgap_tail_profile_test.sh" "$build"
java -ea -Xmx512m -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector \
 -cp "$build/runtime/current" test.prok.ProteinSupportFateTest \
 > "$build/test.prok.ProteinSupportFateTest.stdout" 2> "$build/test.prok.ProteinSupportFateTest.stderr"
cat "$build/test.prok.ProteinSupportFateTest.stderr"
for suite in test.prok.ProteinStopEvaluationTest test.prot.HbmTailChangeAuditTest test.prot.HbmTailLossAuditTest test.prot.HbmTailRobustnessTest test.prot.ReferenceConcordanceReportTest; do
 fixture_args=()
 if [[ $suite == test.prok.ProteinStopEvaluationTest ]]; then mkdir "$build/fate_stop_fixtures"; fixture_args=("$build/fate_stop_fixtures"); fi
 java -ea -Xmx512m -XX:ActiveProcessorCount=1 --add-modules jdk.incubator.vector \
  -cp "$build/runtime/current" "$suite" "${fixture_args[@]}" \
  > "$build/$suite.stdout" 2> "$build/$suite.stderr"
 cat "$build/$suite.stderr"
done
root=$build/fate_cli_cds
[[ ! -e $root ]]; mkdir "$root"
printf '>fixture\nATGAAATAA\n' > "$root/input.fa"
cp "$root/input.fa" "$root/original.fa"
if bash "$build/runtime/callgenes.sh" in="$root/input.fa" out="$root/off.gff" \
 annotation=off hbmdiag="$root/off.tsv" cdsonly=t t=1 tax=f -Xmx256m \
 > "$root/off.stdout" 2> "$root/off.stderr"; then
 echo 'hbmdiag unexpectedly accepted annotation=off' >&2; exit 1
fi
grep -q 'hbmdiag requires annotation=core or annotation=full' "$root/off.stderr"
[[ ! -e $root/off.gff && ! -e $root/off.tsv ]]
if bash "$build/runtime/callgenes.sh" in="$root/input.fa" out="$root/alias.gff" \
 annotation=core hbmdiag="$root/input.fa" overwrite=t cdsonly=t t=1 tax=f -Xmx256m \
 > "$root/alias.stdout" 2> "$root/alias.stderr"; then
 echo 'hbmdiag unexpectedly accepted an input alias' >&2; exit 1
fi
grep -Fq "File '$root/input.fa' was specified multiple times." "$root/alias.stderr"
[[ ! -e $root/alias.gff ]]
cmp "$root/input.fa" "$root/original.fa"
printf 'PROTEIN_FATE_CLI_PASS off_rejected=true input_alias_rejected=true outputs_unopened=true\n' > "$root/PASS"
cat "$root/PASS"

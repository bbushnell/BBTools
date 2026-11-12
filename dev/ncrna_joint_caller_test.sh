#!/bin/bash
# Bounded whole-caller equivalence for both complete grouped families.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: ncrna_joint_caller_test.sh in=cases.tsv reference=registered.tsv sixsmap=rules.tsv cobalaminmap=rules.tsv out=fresh_directory
Case TSV header: family, input_id, fasta, expected_TP (tab-separated).
Requires two qualified original inputs per group, spanning both 6S members and
all three cobalamin members. Actual member calls and work must be nonempty.
Runs fresh ordinary callers, then fixed/adaptive equivalence, a changed complete
group policy, reuse/cache off/on, saved-GFF regrading and summed logical work.
Includes legacy three-column baseline controls. No optimization search is run.
Parent heap defaults to 512 MiB; sequential caller children use 4096 MiB.
Optional members=map.tsv loci=inventory.tsv adds one coverage-enabled sweep per
group, validates complete member accounting, and compares every fixed-tuple GFF
with the ordinary coverage-disabled control. Both metadata files are required.

policycontrols=directory selects the smaller ordinary-policy integration test
instead of the cache/reuse matrix. Supply six original inputs: two each with
family labels 6S, cobalamin and sam1, plus sam1map=rules.tsv. The directory holds
6S.tsv, cobalamin.tsv and sam1.tsv, each containing exactly shipped and selected
complete policy tuples in native nine-column format. It compares no-table and
explicit-shipped ordinary calls, then ordinary versus sweep GFFs and target
logical work for both tuples. All competing callers remain enabled.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" test.prok.NcrnaJointCallerTest "$@"

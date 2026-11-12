#!/bin/bash
# Describe saved call/reference geometry without changing grades.
set -eo pipefail
usage(){
	cat <<'HELP'
Usage: dev/summarizencrnacalloverlap.sh in=grade.raw.tsv family=glycine out=new_prefix

Reads the native GradeNcrnaReference raw report, one original input at a time.
Reports each original FP's largest same-family/contig/strand scored-reference
overlap, mask overlap, and nearest originally matched call with signed gap.
Coordinates are one-based inclusive; gap zero means immediately adjacent,
positive means intervening bases, negative means overlapping calls.

The half-call overlap classification is descriptive only. It changes no TP,
FP, FN, reference coverage or precision number. Exact duplicates remain visible.
Witness multiplicity counts all evaluated calls at least half covered by the
selected reference, including masked calls. Several witnesses may overlap;
this is not a two-aptamer assignment or proof of independent support.

out=prefix writes .calls.tsv, .summary.tsv and .distances.tsv.
The distance table has cumulative 0/10/25/50/100/200/500-base counts for original
FPs with zero scored-reference overlap and a non-overlapping original TP nearby.
It does not select a biological adjacency threshold. Inputs without any selected
family rows are absent from this raw format and are not inferred from filenames.
Model names are retained; numeric verifier identity is unavailable (NA).
No caller, alignment, sequence read, admission or biological adjudication occurs.

Fields preserve native report UTF-8/percent escaping. Family must match the report
field exactly. Rows for each input/family must be contiguous. Existing outputs
are refused; overwrite and append are unsupported. One processing thread.
HELP
}
if [[ $# == 0 || "$1" == -h || "$1" == --help ]]; then
	usage
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	RESOLVED_SCRIPT="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$RESOLVED_SCRIPT"
done
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.SummarizeNcrnaCallOverlap "$@"

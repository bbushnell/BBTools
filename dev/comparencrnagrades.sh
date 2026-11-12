#!/bin/bash
# Fixed-reference caller comparisons using completed native grading reports.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 9, 2026

Usage: dev/comparencrnagrades.sh in1=baseline_grade in2=candidate_grade \
       familymap1=baseline_families.tsv familymap2=candidate_families.tsv out=new_prefix

Compares saved GradeNcrnaReference results without calling or matching genes.
The inputs are grade prefixes, each with .raw.tsv, .inputs.tsv and .summary.tsv.
Only caller settings may differ. Input IDs/FASTA/reference paths, reference locus
identities/geometry/models/dispositions, overlap settings, qualifications and
declared family maps must agree. Reference-swap comparisons are unsupported.
The original family maps must be supplied; declarations do not prove historical
producer provenance. Original logs/configs remain the recipe owner's evidence.
Overlap is compared as reported (six decimal places); original configurations
must establish exact settings when differences below that precision matter.

in1=         Baseline grade prefix.
in2=         Candidate grade prefix.
familymap1=  Original baseline grading family map.
familymap2=  Original candidate grading family map; content must equal familymap1.
out=         Fresh comparison prefix. Overwrite and append are unsupported.
t=1          Native I/O threads; comparison is single-threaded.
-Xmx1g       Default heap, with 256 MB initial allocation.

Outputs:
.loci.tsv    Every input/family/reference locus, both outcomes and signed offsets.
             Endpoint-change flags apply only to TPs retained in both runs.
.counts.tsv  Original per-family summaries tagged baseline/candidate; preserves
             TP/FP/FN, masks, duplicates, exclusions, precision and recall.
.summary.tsv Retained/lost/gained TPs, shared FNs, masks and retained-TP endpoint
             changes. Detection losses do not also count as lost exact endpoints.
             Written after detailed outputs close successfully.

Counts and query/reference witnesses are checked before outputs open. This detects
incomplete/mismatched reports; it does not repeat biological validation or certify
that the original caller configs accurately describe the saved calls. Encoded text
fields retain the native grader's %HH escaping and are compared without decoding.
HELP
}
if [ $# -eq 0 ] || [ "$1" = -h ] || [ "$1" = --help ]; then
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
parseJavaArgs --xmx=1g --xms=256m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.CompareNcrnaGrades "$@"

#!/bin/bash
# Diagnostic comparison only; does not evaluate biological admission predicates.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 8, 2026

Usage: dev/comparencrnascenarios.sh in=before.tsv in2=after.tsv calls=calls.tsv \
       contigs=contigs.tsv familymap=families.tsv out=new_directory

Compares two explicitly supplied hypothetical dispositions against the same saved
CallGenes calls. No searches, biological predicates, cutoff selection or approval
claims. No reference GFFs or reference manifest are emitted. Changes in TP/FP/FN
are consequences of different assumed references, NOT measured accuracy gains.

Inputs have tab-separated headers (shown with spaces below):
in=/in2= decision_id input_id literal_target start stop strand reference_family family_state disposition rule_id model reason decision witnesses
         The thirteen disposition fields follow dev/writencrnareference.sh's schema,
         with an explicit stable decision_id added first. Key=(input_id,decision_id)
         must be unique within each scenario. Shared keys declare correspondence;
         they do not independently prove biological identity. Changed coordinates
         or families stay paired. Different keys remain additions/removals; no
         overlap-based split/merge correspondence is inferred.
         Text uses evidence TSV escaping (%25, %09, %0D, %0A); status tokens are
         literal. Unknown model/coordinates use a dot. family_state is RESOLVED or
         UNRESOLVED; disposition is scored/dont_care/unresolved/removed/context.
         Rule, decision and witness pointers must be nonempty. Non-scored rows
         need a reason. A rule identifier is provenance, not an APPROVED claim.
calls=   input_id fasta gff
         Literal UTF-8 fields and absolute original FASTA/call-GFF paths. Input
         IDs and paths must be unique. Include empty-call inputs explicitly.
         This one roster binds BOTH scenarios; scenario rows cannot add inputs.
contigs= input_id fasta literal_target observed_length
         Escaped saved observations, exact roster FASTA paths, unique targets.
         Every declared input needs observations. Generic calls and scored/masked
         scenario rows must fit the corresponding observed contig bounds.
familymap= family model_prefix display_label include
         Same explicit literal-UTF-8 CallGenes projection as dev/gradencrnareference.sh.
         include is true/t or false/f. Scored/masked families must be included.

Only resolved scored/masked dispositions participate in hypothetical matching.
Unresolved/removed/context rows remain in the accounting without scoring authority.
Duplicate/contradictory resolved physical dispositions fail. Matching is one-to-one;
scored eligibility takes precedence over same-family masks and duplicate handling.
Reference order is normalized before matching. Every original call is reported,
including excluded and duplicate calls. Structural RNA feature counts are separate.

Outputs in a fresh directory:
sources.tsv         Original diagnostic input paths (not a reference manifest).
dispositions.tsv.gz Every original row with its scenario and source line.
changes.tsv         Declared-key pairs/additions/removals, original source lines,
                    both intervals/families/models/statuses and change flags.
calls.tsv.gz        Both hypothetical outcomes, decision/physical-ID witnesses,
                    source lines, duplicate/eligibility flags and changed flag.
inputs.tsv          Exact original FASTA/call bindings and feature/decision counts.
counts.tsv          Before/after disposition counts by family and model.
summary.tsv         Before/after family TP/FP/FN, scored/masked denominators,
                    masked/duplicate/excluded calls and signed after-minus-before
                    TP/FP/FN/masked-call deltas and effective overlap setting.
                    Written last; require exit zero.

Every report is labeled HYPOTHETICAL_REFERENCE_COMPARISON_NOT_BIOLOGICAL_ACCURACY.
Declared input bindings and correspondence are not independent producer evidence.
Unresolved physical coordinates remain in provenance; no interval is repaired.
Defaults: overlap=0.5 reciprocal scored overlap, t=1, -Xmx4g, ow=f.
Unicode filenames require a supporting native filesystem locale. Existing output
paths are refused. No model/corpus/FASTA payload is read and no calling is run.
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
parseJavaArgs --xmx=4g --xms=256m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.CompareNcrnaScenarios "$@"

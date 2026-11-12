#!/bin/bash
# Serializes explicit reviewed disposition claims; does not decide biological truth.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 9, 2026

Usage: dev/writencrnareference.sh in=dispositions.tsv inputs=inputs.tsv \
       contigs=evidence.contigs.tsv.gz rules=reviewed_rules.tsv out=new_directory
       dev/writencrnareference.sh mode=rebase in=reference.tsv out=new_gff_directory \
       evidencefrom=/old/annotation_witnesses.tsv.gz \
       evidenceto=/keeper/annotation_witnesses.tsv.gz

mode=write (default) serializes dispositions as described below.
mode=rebase relocates an existing reference GFF view. The input manifest has the
literal header input_id, fasta, gff (tab-separated). Input IDs and original FASTA
paths remain unchanged; output GFF paths point into the fresh directory.
The destination witness must already exist outside that fresh output directory.
The owner copies/verifies witnesses and other reconstruction inputs separately.

Rebase changes only each feature's reference_evidence value, requiring an exact
match to the escaped evidencefrom path. It preserves other line content, including
unknown attributes, comments, numbers, masks and empty views; line endings become
LF. For canonical LF output, decompressed bytes differ only at that value.
Compression bytes need not match. Source files are never overwritten.
File/input directives and explicit row input/status/evidence fields are checked;
this is not a biological, full-schema or witness-copy validation. Raw source GFFs
without the required reference fields fail. Inputs/rules/contigs options apply only
to mode=write; evidencefrom/evidenceto apply only to mode=rebase.
Rebase produces input_N.gff.gz, rebase.tsv, and reference.tsv published last.
An incomplete directory without reference.tsv is not a completed reference view.

Mechanically writes one reference GFF per declared input, including empty files.
It does not select thresholds, rank models, merge loci, infer families or repair
coordinates. The input disposition table must already contain those decisions.
A separate rule table must record APPROVED and an authority pointer for every
rule. These are required claims, not independent verification of human approval
or biological truth. A scored word on a locus does not bypass the rule table.

All schemas below have a literal tab-separated header, shown here with spaces:
inputs=  input_id fasta
         Literal UTF-8 fields; absolute original FASTA paths. IDs and exact paths
         must be unique and contain no tab/CR/LF. FASTA payloads are not read.
contigs= input_id fasta literal_target observed_length
         Existing saved contig observations. Paths must match the roster exactly;
         targets are unique within each input. Every input needs observations.
rules=   rule_id family disposition approval authority
         approval must be APPROVED; authority must be nonempty and not a dot.
         disposition is scored, dont_care, unresolved, removed or context.
         This is a serialization-rule table, not the provisional model-cutoff draft.
in=      input_id literal_target start stop strand reference_family family_state disposition rule_id model reason decision witnesses
         family_state is RESOLVED or UNRESOLVED. Family and disposition must match
         the explicit approved rule. Non-scored rows need a reason; every row needs
         decision and witness pointers. Dot denotes absent model/unknown coordinates.

Text fields in contigs/rules/dispositions use the evidence TSV escaping for percent,
tab, CR and LF (%25, %09, %0D, %0A), preserving UTF-8. Structural status tokens are
literal. Source row contents are preserved in the companion, with LF line endings.

Only RESOLVED scored/dont_care rows with valid observed contig bounds and explicit
strand enter GFFs. A mask requires a named resolved family and a reason. Other
states remain companion-only, including unassessable original coordinates. For
resolved valid physical identities, duplicate or contradictory dispositions fail;
no winner is selected. Original whole intervals remain unchanged.

out=     Fresh directory without tab/CR/LF; percent and Unicode are supported.
         Unicode filenames require a supporting filesystem locale (for example
         LC_ALL=C.utf8 on Linux); unrepresentable names are rejected before creation.
         Produces input_N.gff.gz files, rules.tsv,
         dispositions.tsv and summary.tsv. The grader-compatible reference.tsv
         is atomically published LAST, after all files close successfully. Its
         input_id/fasta/gff fields are literal UTF-8. A partial directory without
         reference.tsv is not a completed reference set. Existing paths are refused.

Defaults: t=1, -Xmx4g, ow=f. Computation is metadata-only and processes one stream;
BBTools compression flags remain available. No caller or model search runs.
An empty emitted view is not proof of biological gene absence. This command does
not implement admission predicates or approve applying them to a real dataset.
Diagnostic scenarios must not be presented as an authoritative biological reference.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.WriteNcrnaReference "$@"

#!/bin/bash
# Grades declared per-input calls against explicit scored/masked reference GFFs.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 9, 2026

Usage: dev/gradencrnareference.sh in=queries.tsv reference=references.tsv \
       familymap=families.tsv out=new_prefix

Grades CallGenes queries one input at a time using explicit reference_status.
It does not select biological thresholds, call genes, admit reference loci, or
verify that producer input claims are true. Do not relabel a combined saved GFF
as one input: retain its original input mapping and split provenance first.

Required:
in=          Query manifest with literal TSV header: input_id fasta gff
reference=   Reference manifest with the same header and, by default, the same inputs.
             Both use literal UTF-8 fields and absolute FASTA/GFF paths. FASTA
             paths must agree exactly; sequence payloads are not opened. Duplicate
             IDs, FASTA paths or reused GFF paths fail. An input with zero calls
             still needs an explicit query row pointing to an empty GFF file.
familymap=   TSV header: family model_prefix display_label include
             include is true/t or false/f. Every generic RNA model needs exactly
             one explicit mapping; exclusions remain counted in the reports.
             Several prefixes can project to one biological scoring family.
             Present CallGenes family labels must corroborate the model mapping.
out=         Fresh output prefix. Existing outputs are refused.

restrictreference=f Default strict mode requires equal complete input populations
             and an included projection for every reference family. Set true to
             select the complete query input population and included families
             from a larger reference. Every query input still needs an exact
             reference ID and matching FASTA; missing/duplicate bindings fail.
             Off-panel GFFs are not read. Empty query files keep eligible FNs.
             Selection follows the declared family map, never observed calls.
             The family map/config and complete query manifest must describe the
             original baseline; this flag does not independently prove that claim.

Options:
overlap=0.5  Minimum reciprocal overlap for scored matches; must be in (0,1].
             Masks use any positive same-input/family/contig/strand overlap only
             when the call has no eligible scored reference. Exact duplicate
             evaluated calls never consume additional scored references.
t=1          Native I/O thread setting; grading processes one input at a time.
-Xmx4g       Heap; launch defaults are 4 GB maximum and 256 MB initial.
ow=f         Fresh files required. Overwrite and append are unsupported.

Optional cobalamin subtype reporting (all three arguments required together):
subtypes=    TSV header: subtype model_prefix reference_model. Unique names,
             prefixes and reference models; covers every included cobalamin rule.
witnesses=   Registered annotation_witnesses.tsv.gz with parent/event bindings.
events=      Its event_provenance.tsv.gz with retained tier decisions.
             The project must authenticate these consumed tables against its
             pinned keeper manifest. This tool checks bindings, not authorization.
             Protected prior subtype assignments remain fixed. Weak and bridge
             observations do not confer support; marginal support stays unscored.
             Scored parents without justified scored subtype support fail loudly.

Both subtype views are computed in ONE invocation from the same loaded loci:
MODEL_SUPPORT_RECOVERY grades each supported subtype independently. Its FN means
model-support recovery, not a missed gene, not wrong-subtype proof.
EXCLUSIVE_SUBTYPE_ACCURACY masks multi-supported parents in their supported
subtypes before matching. Unsupported third subtypes receive no ambiguity mask.
This is a reference-support restriction, not proof of unique biological subtype.
Existing union masks, scored-before-mask behavior and duplicate accounting remain.
Subtype totals are nonadditive. The ordinary union grade is unchanged.

Reference files must use the NcrnaReferenceGff version/input directives and explicit
scored or dont_care rows, with reasons for masks. Unknown status and duplicate locus
IDs fail. Queries use CallGenes' RNA features and comma model:/family: attributes;
CDS, tRNA and rRNA rows are counted separately. This is not a generic GFF grader.

Create reference= with a public serializer, rather than passing a source GFF:
  dev/writencrnareference.sh in=dispositions.tsv inputs=inputs.tsv \
    contigs=contigs.tsv rules=reviewed_rules.tsv out=reference_directory
Then use reference=reference_directory/reference.tsv. The writer's help gives
the complete input headers and approval/disposition fields. Alternatively,
dev/buildncrnareference.sh mode=tiers produces reference.tsv from completed bound
CM evidence and explicitly selected tier rules; its help describes those inputs.
Both producers write each reference GFF with these first two lines:
  ##gff-version 3
  ##reference-input-id ORIGINAL_INPUT_ID
The input ID must agree with the manifest. Feature attributes include input_id,
reference_family, reference_status and reference_evidence; masks also require
dont_care_reason. The serializer constructs the coordinate-qualified ID, which
the grader checks; a free-form locus ID cannot replace it.

Outputs:
.raw.tsv     Every reference and generic-query row, outcome and matching/mask
             witnesses. Excluded queries remain visible; their duplicate status
             is unassessed. d5/d3 apply to matched queries, in transcript direction.
.inputs.tsv  Exact declared paths and included/excluded/non-target feature counts.
.summary.tsv Per-family TP/FP/FN, scored/masked counts, evaluated duplicates,
             exclusions, exact endpoint counts, precision and recall. Zero
             denominators are NA. No genomic false-positive rate is estimated.
             Written only after all inputs and detailed outputs complete.
.selection.tsv Only with restrictreference=t: retained/omitted INPUT counts plus
             kept/omitted scored/masked ROW counts by family among inspected
             on-panel GFFs. Off-panel row totals are unmeasured, not zero.
             Unmapped reference families and explicitly excluded query families
             are distinguished. Original global reference files remain unchanged.
.subtype_raw.tsv Optional: complete separately matched, view-tagged outcomes.
.subtype_support.tsv Optional: each parent once, original bounds/status, scored
             and marginal subtype support, protected authority and ambiguity flag.
.subtype_summary.tsv Optional: both views' detection, masking and endpoint counts.
.subtype_counts.tsv Optional: distinct parent counts, including multi-supported.
             Choose a presentation by selecting a stored view, never by filtering
             recovery TP/FN labels to derive exclusive assignments.

Report text escapes percent/tab/CR/LF as %HH and otherwise preserves UTF-8.
Reference masks remove uncertainty from scoring; they can hide actual false calls.
A successful grade validates this computation, not the biological reference.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.GradeNcrnaReference "$@"

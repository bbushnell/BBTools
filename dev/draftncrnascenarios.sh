#!/bin/bash
# Explicit hypothetical rules over saved evidence; never a biological reference.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden
Last modified November 8, 2026

Usage: dev/draftncrnascenarios.sh in=evidence_prefix census=genomes.tsv \
       memberships=memberships.tsv scope=scope.tsv assets=assets.tsv queries=models.tsv \
       policy=draft.tsv priors=bound_priors.tsv.gz links=source_links.tsv.gz out=new_directory

Constructs two hypothetical dispositions using saved dev/buildncrnareference.sh evidence.
The first applies structural checks at each model's GA. The second additionally
requires the explicit draft margin for otherwise eligible models. Nothing is
adopted as biological truth. No scans, caller runs, FASTA reads, reference GFFs or
reference manifests are produced. Counts are effects of assumptions, not accuracy.

in= Saved GA-filtered evidence: summary, coverage, sources, hits, overlaps, contigs.
    Every native event must pass its captured model's GA, including equality, and
    its saved GA flag must agree. Below-GA rows are rejected, not silently dropped.
    This interface does not admit an unfiltered sensitive scan. Prior positives
    remain independent of current score or absence from the native scan.
census=/memberships=/scope=/assets=/queries= The same declared catalog/model tables
    used for that evidence. No model or genome payload is opened.
links= The complete saved strata.source_links.tsv.gz table. Direct same-strand
    RFAM claims participate for every row kind, including regulatory/tagged genes;
    pseudo/partial/exception attributes affect routing only on RNA rows themselves.
    All original links remain
    in the original companion table, joined by input/native-table/native-line.
priors= Complete bound prior_decisions.tsv.gz from mode=prior, including original
    source/core/whole coordinates, classification and authority. Exact accepted
    cobalamin cores and accepted novel FMN whole intervals are preserved, even
    without a new hit. Other prior rows gain no inherited masking exemption.
policy= Explicit twenty-three-column draft table (literal tab-separated header):
    accession family release model_sha80 GA_bits model_length prior_anchor_N
    prior_score_N distribution_basis anchor_N score_N full_score_N full_margin_p05
    full_margin_median full_margin_p95 coverage_p05 coverage_median
    illustrative_p95_margin_bits new_candidate_disposition qualification status
    anchor_tier distinct_input_N
    status must be DRAFT_NOT_APPLIED and captured model identity must agree.
    BOUNDARY_UNSELECTED compares GA against the declared nonnegative p95 margin.
    NO_AUTOMATIC_ADDITION keeps new candidates unresolved.
    CONTEXT_PENDING_BRIAN keeps model evidence as context, never mature tRNA truth.
    MASK_PENDING_SPARSE_ANCHOR_REVIEW, MASK_PENDING_STRUCTURE_REVIEW and
    MASK_PENDING_SUBTYPE_REVIEW remain fixed masks when family is unambiguous.
    The latter actions require NA in illustrative_p95_margin_bits. Queried models
    without a policy remain context and still participate in competition checks.

Rules fixed in both scenarios:
- All valid, GA-passing queried-model same-strand positive overlaps participate. Exact same-model
  interval duplicates collapse; conflicting duplicate scores/spans stay unresolved.
- Different-family competition or direct unknown/other-family RNA claims remain
  unresolved, never one mask per competing family. A preserved positive prior
  takes precedence: overlapping new candidates supply provenance, no wider mask.
- Same-family competing intervals, direct RNA source warnings, partial model span,
  native truncation and physical contig edges block auto-scoring.
  Their family assignment is an explicit diagnostic assumption, not adjudication.
- Direct RNA pseudo/partial/exception/unrecognized-boolean/incompatible-ID warnings
  and different-accession same-family claims block auto-scoring. Context/component
  attributes do not transfer to an RNA event.
- Masks collapse only at identical input/family/contig/strand/bounds. All original
  model witnesses remain; a multi-model physical mask uses model dot. No transitive
  merging or selected highest-scoring duplicate is performed.

Outputs (new directory, no overwrite):
ga.tsv.gz / margin.tsv.gz  Explicit comparator-schema hypothetical dispositions.
    Stable input-local event IDs declare correspondence, not biological identity.
witnesses.tsv.gz  Every original event and its emitted decision, model-interval
    representative, reasons, flags and required source-link companion join.
    Flag bits:1 cross/unknown RNA claim,2 RNA warning,4 conflicting duplicate,
    8 same-family competitor,16 other-family competitor,32 preserved-prior overlap.
prior_witnesses.tsv.gz Every unchanged prior row with its generated decision ID.
counts.tsv       Unique model-interval dispositions per model/family/scenario.
family_counts.tsv Emitted rows per family/scenario, including priors and counting
                 coincident masks once. Unresolved rows need not be distinct loci.
summary.tsv      Raw events, duplicate events, prior/policy counts, observed-input
                 scope and emitted disposition rows; written last, require exit0.
sources.tsv      Original input paths and the hypothetical-result qualification.

Raw events, model intervals and emitted physical masks are different counts.
Unresolved/context rows may describe the same coordinates without claiming genes.
All source-link/native evidence files remain required companions; do not discard
those original tables after construction. Subset these global scenario tables to
EACH saved baseline's exact inputs and historically enabled/supported families
before comparison. Complete empty-query inputs still contribute enabled-family FNs.
The tool does not infer baseline scope, choose a cutoff, or certify a family.
Defaults: one analysis worker, -Xmx4g, fresh outputs, assertions enabled.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.DraftNcrnaScenarios "$@"

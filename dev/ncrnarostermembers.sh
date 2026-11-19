#!/bin/bash
# Native construction members from a complete reconciled ncRNA roster.
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
	cat <<'HELP'
Usage: dev/ncrnarostermembers.sh in=reconciliation_result cohorts=accepted_cohorts.tsv out=fresh_directory expectedloci=N t=16 -Xmx4g -eoom
Reads merged_loci.tsv.gz, source_frames.tsv.gz, locus_components.tsv.gz,
loci.tsv.gz, events.tsv.gz and summary.tsv from the accepted reconciliation.
cohorts= is the seven-column input manifest, not its five-column result summary.
Construction uses native cmsearch bounds/strand only, with no old annotation
endpoints or old library sequences. Reads the recorded canonical source path;
a supplied symlink is provenance and cannot redirect construction.
Run mode=metadata first: reports eligibility.tsv, eligibility_calls.tsv.gz,
eligibility_parameters.tsv and missing_model_positions.tsv without reading FASTA.
Construction excludes Infernal-truncated hits. coveragefilter=f (default)
retains incomplete-model-coverage hits with their own cmsearch bounds. Coverage
stays diagnostic: model start>modeltolerance OR model end<length-modeltolerance
(default5). Contig-edge contact alone does not exclude. Explicit coveragefilter=t
enables the additional coverage exclusion and maxcoveragefraction=0.25 stop;
that fraction stop is inactive when coveragefilter=f. Member headers retain
native cmfrom/cmto/cmlen and incomplete_model_coverage for cluster diagnostics.
All original calls, including excluded partials, remain in recall.
The same FASTA pass records model-oriented available flanks and immediately
adjacent N-run lengths. genomic_context.tsv reports contig-end reach (missing
model positions exceed available genomic bases), adjacent gaps, internal calls,
other boundary states and unresolved calls. This is an approximate reach test,
not a new alignment. assembly_loss.tsv counts groups left with no original calls
after each proposed removal; bit_scores.tsv reports retained native-score quartiles.
Extracts inclusive bounds in the original strand orientation with native FASTA
streaming. Each positive source frame is read once. Construction duplicates
require the same accepted assembly group and equal complete case-normalized
DNA IUPAC arrays; different groups, lengths and nucleotide symbols remain distinct.
Identical sequences across accepted groups are counted and kept separately.
Writes construction members.fa.gz, merged_members.tsv.gz, original_calls.tsv.gz,
summary.tsv and parameters.tsv. Every original called coordinate group remains
linked through its overlap component and construction representative, or carries
an explicit unresolved extraction reason. Missing/ambiguous contigs do not vanish.
No clustering, model construction, holdouts or cross-group alias resolution.
Default mode=full forbids sampling. Explicit mode=canary maxsources=N retains
the full input ledger but labels skipped frames; that output is operational only.
Review unresolved counts before using a full corpus for model construction.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [[ -h $SCRIPT ]];do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=4g --xms=256m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.NcrnaRosterMembers "$@"

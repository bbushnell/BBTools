#!/bin/bash
# Complete original-gene baseline grading and all-data selected-window GC fit.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/summarize_ncrna_roster_baseline.sh in=capture_prefix population=prepared_prefix assignments=new_builder.tsv out=fresh_prefix expectedgenes=N expectedframes=N expectedmembers=N t=1 -Xmx2g -eoom
For disjoint captures, replace in= with captures=manifest.tsv. Its header is:
capture_id<TAB>capture_prefix<TAB>reference_manifest
Use unique simple capture names and absolute prefix/reference paths. Each native
slot list must exactly match its ordered shard reference and the global ledger.
Duplicates, missing inputs (even zero-call inputs) and cross-shard substitutions
fail before complete results. Captured/completed event IDs then have the form
capture_id:iNeM or capture_id:iNsM; their suffix is the native local event ID.
One global saved-data pass performs matching and GC fitting across all shards.
Reads saved native inputs, completed window results and input-bound final calls.
Every original gene stays in recall. Missing input, no window, unaccepted and
ambiguous-best proposals remain explicit. Membership uses the accepted own span,
not overlap of its larger window, then highest identity and deterministic ties.
An equally best multi-gene proposal leaves membership explicitly unresolved.
For an explicitly requested partial-input sample, add samplemanifest=reference.tsv
with input_id, fasta and gff columns and the sample expectedgenes/expectedframes.
All original genes in those exact frame/path bindings remain in the denominator;
the output scope and success marker explicitly say SAMPLE, not full population.
This mode is unavailable for acceptance-only parity controls.
Final detection is separate maximum-cardinality one-to-one reciprocal1/2 matching.
Contig first tokens and ordinals must bind one-to-one across both saved events
and final calls within each source; conflicting bindings fail loudly. Previously
extractor-unresolved genes remain explicit recall misses and cannot gain GC
membership from a matching name. Construction-only exclusions remain in recall.
Unmatched final calls are not adjudicated false positives. Counts are native event
work summed once, independent of overlapping gene associations. Model GC uses the
selected actual window, ACGT-only denominator, population SD and explicit absent/
zero-SD states. No DNA, alignment, caller or rank replay. All data, no holdouts.
Outputs .genes.tsv.gz, .gene_windows.tsv.gz, .calls.tsv.gz, .model_gc.tsv and
.summary.tsv. Per-gene overlapping/selected event work and all gene/window
associations retain actual counted work; associated costs are nonadditive across
genes. Global cost counts each event once. Completed event IDs use input slot and
within-input completed-row ordinal; captured IDs are the original window FASTA IDs.
Associations retain reference-resolution state separately from geometric overlap.
For real saved-window replay, add acceptanceonly=t identity=0.60|0.62|0.64.
Use its replay prefixes as in= or captures= entries, retaining the same original
population and construction assignments. At64%, parity=original.genes.tsv.gz is
required and checks the exact accepted original-gene set, including every miss.
Acceptance-only mode reads no final calls, fits no GC model, and writes no
.calls.tsv.gz or .model_gc.tsv. It writes .genes.tsv.gz, .gene_windows.tsv.gz,
.cohorts.tsv and .summary.tsv with actual alignment work, coverage/truncation
subsets, and unadjudicated accepted events without roster overlap. Per-window,
per-million-saved-window-base and per-million-recorded-input-base denominators
are explicit; overlapping cohort counts and per-gene work are nonadditive.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ];do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.SummarizeNcrnaRosterBaseline "$@"

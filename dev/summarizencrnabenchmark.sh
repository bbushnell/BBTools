#!/bin/bash
# Read-only timing and exact-parent comparison of native ncRNA grade reports.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell, Raiden and Ganyu

Usage: dev/summarizencrnabenchmark.sh in=times.tsv gradeinputs=grade_inputs.tsv \
       grades=grade_summary.tsv raw=grade_raw.tsv out=new_prefix \
       policies=shipped,selected workers=1,16 repeats=3 baseline=shipped

Timing without regrading: dev/summarizencrnabenchmark.sh timingonly=t in=times.tsv \
       out=new_prefix policies=shipped,selected workers=1,16 repeats=3

Summarizes completed ordinary-caller benchmark reports. Reads no sequences and
does no calling, alignment, regrading or reference admission. Requires the exact
declared policy x worker-count x repeat population, including zero-query inputs.

Inputs:
timingonly=f Explicit timing-only mode rejects gradeinputs/grades/raw/counts/taxonomy.
             It also rejects subtype reports.
             Writes only cells/timing and optional profile; accuracy and detection
             stability are NOT evaluated. No grade/comparison tables are produced.
in=          Header: policy threads repeat user_seconds system_seconds wall_seconds
gradeinputs= Native GradeNcrnaReference .inputs.tsv with policy/threads/repeat
grades=      Native .summary.tsv with the same three leading fields
raw=         Native .raw.tsv with the same three leading fields
counts=      Optional binder original_counts.tsv with those leading fields
taxonomy=    Optional literal TSV: input_id fasta domain phylum corpus_inputs
             category taxonomy_status. Exact original IDs/FASTA paths must bind
             every selected input once, including zero-query inputs. A larger
             metadata population is allowed; off-panel rows are not validated.
             Labels and corpus counts are declarations, not inferred taxonomy or
             a new rarity threshold. Unknown labels stay explicit. Fields are
             literal UTF-8; native grade tables retain their percent escaping.
subtyperaw=   Optional native cobalamin .subtype_raw.tsv with policy/threads/repeat
subtypegrades= Its .subtype_summary.tsv with the same three leading fields.
             Supply both or neither. Every declared cell needs all three existing
             cobalamin subtypes in both native views, including zero counts.
             Subtype totals, input populations, parent geometry/masks and actual
             query-line/model witnesses are checked against the ordinary reports.
             A net TP decline is not a lost-set size. No subtype is inferred from
             reference coordinates or from the model winning a different match.
profile=     Optional native NcrnaCpuProfile table (requires profiletime=)
profiletime= That SAME profile invocation's single user/system/wall time row
             The profile has no policy label; the recipe binds it to its policy.
policies=    Required comma-separated literal policy names, no empty/duplicate names
workers=1,16 Declared positive worker counts; these are data labels, not tool threads
repeats=3    Declared positive number of fresh-JVM observations per cell
baseline=shipped  One of the declared policies

Native grade fields retain their original percent escaping and UTF-8 bytes.
Compared cells must share inputs, exact FASTA/reference paths, scored/masked
reference-parent keys and overlap rule. Counts and endpoint totals must reconcile
with raw references and queries. This checks declared reports, not provenance or
genome contents. Do not relabel a combined GFF as a graded input.
The shipped one-worker baseline needs at least 10 CPU seconds, following the
cluster timing rule. Invalid numeric fields fail through BBTools assertions;
launchers keep assertions enabled.

Outputs:
.cells.tsv        Per-repeat seconds and optional original feature counts
.timing.tsv       User/system/CPU/wall medians and ranges; matched-worker savings
                  CPU is user+system PER REPLICATE before the median.
                  Timing-only rows state that caller outputs were not compared.
.grades.tsv       Reconciled native counts and ratios, not new scoring
.comparisons.tsv  Exact lost/gained parents and first/last detected inputs per family
.changes.tsv      Exact input/family/reference keys for each lost/gained detection
                  Policy vs baseline uses matching workers/repeat. Additional
                  comparisons check repeats and worker counts of the same policy.
.lineages.tsv     With taxonomy=: per declared class, selected-input denominator,
                  TP/lost/gained references and detected/first/last input counts.
                  Class totals reconcile with every family comparison, including
                  classes with no detections. Labels never change a native grade.
.input_changes.tsv With taxonomy=: exact first-detection/last-loss input IDs and
                  their labels, declared corpus count and before/after TP counts.
                  Replacing a copy within an already detected input is neither.
.profile.tsv      Optional disjoint worker CPU over its own process CPU denominator
                  Shared CDS/RNA path is excluded from RNA-only totals. Resource
                  loading and other-thread work remain unallocated. Clock overhead
                  remains included; profile timing is not a performance cell.
.subtype_comparisons.tsv  Optional exact subtype-parent gains/losses and counts by
                  ordinary union outcome/model. Both views remain nonadditive.
.subtype_changes.tsv  Optional changed input/parent keys with coordinates, both
                  subtype outcomes and both ordinary union outcomes/query lines/
                  models. LOST assesses the right union; GAINED assesses the left.
                  OTHER_SUBTYPE_MODEL and SAME_SUBTYPE_REMATCH mean the ordinary
                  union still matches that exact parent with the indicated model.
                  UNION_MISS means its ordinary outcome is FN; UNION_MASKED means
                  DONT_CARE. These are saved matcher outcomes, not biological
                  absence, unique subtype identity, or a new matching decision.
                  All existing reports remain unchanged by these optional inputs.

Fresh outputs only; ow=f and append=f. Exit zero and the final completion report
are required; a parsing or I/O failure can leave partial outputs. This tool does
not prove that reference keys represent distinct biological genes.
Full-analysis completion counts differing ordinary family-union repeat and worker
comparisons once per cell comparison. Subtype-view changes are reported separately
and do not increment those counters; zero there does not certify subtype stability.
Subtype and union views match independently: a subtype TP need not be a union TP
for the same parent. The original query can be assigned to a different parent in
the other view, which is why both outcomes and query witnesses are retained.
Timing-only mode
uses NCRNA_BENCHMARK_TIMING_ONLY_COMPLETE and withholds stability claims.
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
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.SummarizeNcrnaBenchmark "$@"

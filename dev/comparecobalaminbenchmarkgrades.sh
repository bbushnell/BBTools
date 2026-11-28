#!/usr/bin/env bash
# Saved-assignment comparison and overlap-only presence diagnostics.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/comparecobalaminbenchmarkgrades.sh in1=export_before in2=export_after out=fresh_prefix
Optional family=cobalamin selects the exact canonical family used for both exports.
Other families require matching .family.tsv identity files and report one biological
union cohort, including first/last-detection genomes and fixed paired endpoints.
Their version markers begin NCRNA; default cobalamin report bytes are unchanged.
Member-label transitions are not measured for these unions: their summary rows
are omitted and paired label_switches fields are NA. Any present .family.tsv
declaration is checked even with default or explicit family=cobalamin.
The changesonly, subset, replay and mode=report options remain cobalamin-only.
Presentation: dev/comparecobalaminbenchmarkgrades.sh mode=report in=completed_tables \
  in2=original_replay out=fresh_directory t=1
Report mode reads the six small-radius configurations on 1000/300 genomes,
checks table schemas and unique member/end rows, and writes report.html. It
preserves NA, separates absent from overlap-present losses and keeps fixed paired
cohorts. It does not call, grade, fit or select. Endpoint PNGs are a separate plot step.

Inputs are complete ExportCobalaminBenchmarkGrade views from the SAME reference.
Optional changesonly=t writes only .unmatched.tsv, .unmatched_links.tsv,
.lost.tsv, .lost_links.tsv and .change_summary.tsv. Each unmatched call has one
old-outcome disposition, all positive overlap links are retained, and lost scored
references form a separate population. Geometry is observed, ancestry is not.
Historical per-call CM scores/contact/anchor flags are unavailable in these
inputs. This option cannot combine with radius, trace or subset reporting.
Optional subset=main.fa,minor.fa,third.fa reports named-model retention only when
both libraries come from the SAME original generation and candidate is a subset.
Optional radius=6 adds .paired_endpoints.tsv and .paired_loci.tsv for reference
IDs matched in BOTH views. Reports full paired RMS in bases, better/worse/unchanged
absolute errors, better/(better+worse), and secondary RMS on the SAME ends whose
baseline absolute error is within radius. Empty denominators produce NA.
Member cohorts use the baseline emitted label; switches remain explicit.
Also writes .paired_details.tsv on that same paired N: exact and within-three
counts before/after, RMS and mean absolute distance before/after, moved ends,
and better/worse/unchanged counts. Movement compares signed errors; a +3 to -3
change is moved but has unchanged absolute distance. Empty populations have no
detail statistics; their zero N remains in .paired_population.tsv.

Optional trace=saved_ungated.trace.tsv.gz bindings=panel.sources.tsv with radius=
adds .replay_loci.tsv, .replay_endpoints.tsv, .replay_accounting.tsv and
.replay_selection.tsv. Exact source/contig/strand/model/final-coordinate joins
bind every selected trace to its ON query. Zero-margin reconstruction is required.
Replay margins default to 0, .005, .02 and .1. Explicit replaymargins= accepts
an ordered finite list starting at zero. Record when each value was chosen;
using the option does not make a post hoc value predeclared. Scalar second-end scores use the
captured held/moved first-end context. Full paired RMS is primary, OFF-member
cohorts stay fixed and unknown eligibility makes selection unavailable.
This is conditional score replay, not final detection evidence. In the current
experiment ONLY the 1,000 panel selects a margin; old 300 confirms it unchanged.

Writes .summary.tsv, .genomes.tsv, .loci.tsv, .overlaps.tsv and .lengths.tsv.
Existing one-to-one assignments stay unchanged. Lost/gained reference IDs,
emitted-member switches and raw versus scored first/last genome calls stay separate.
All same-input/contig/strand overlaps at changed loci retain model, assignment,
intersection/reference/call lengths and fractions. Presence is not graded recall.
Length tables retain matched and overlapping-unmatched pairs separately; a query
may overlap multiple references. Such pairs are not additional true positives.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.CompareCobalaminBenchmarkGrades "$@"

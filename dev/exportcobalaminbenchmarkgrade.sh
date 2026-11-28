#!/usr/bin/env bash
# Export completed union grades; no caller or matching operation runs.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/exportcobalaminbenchmarkgrade.sh in=completed_grade out=fresh_prefix
Optional family=cobalamin selects the exact canonical family from the grade.
For example, family=RNaseP exports that biological union without inferring subtypes.
Other families use canonical export mode, write one family cohort and a .family.tsv
identity file, and omit cobalamin alignment columns from metrics. Legacy, members,
model-roster and timing/callerlog options remain cobalamin-only.
Generic normalized cohort columns are named family_union; original model fields
stay unchanged, and no emitted subtype is inferred.
For historical CobalaminPanelGrade: legacy=t panel=original_panel.loci.tsv
For a canonical grade produced with inexact=t: inexact=t

Explicit inexact mode reads the recorded .overlap.raw/.overlap.summary reports
and ordinary .inputs roster. It preserves every assigned match, including
endpoints-inexact overlaps, without rematching. Default exports remain strict;
legacy=t cannot be combined with inexact=t. Metrics label the inclusive view.

Writes .members.tsv, .genomes.tsv, .witnesses.tsv, .matches.tsv,
.queries.tsv, .references.tsv and .metrics.tsv. The Java metrics row gives total
calls and scored/masked reference counts separately from masked query calls.
Optional timing=native.time.tsv callerlog=complete.stderr copies actual
user/system/wall seconds and the three cobalamin alignment counters into metrics.
Query/model/coordinates remain available for
overlap-only loss classification and matched versus unmatched length diagnostics.
Emitted members partition UNION outcomes; these are not biological subtype grades.
The nine-column matches contain clan and member views for ncrnaendpointsummary.sh.
Offsets retain gene-sense signs and actual matched reference coordinates.
Zero-call genomes survive. No calling, rematching, DNA read or truth change.

Optional small member summary over one or more DISJOINT canonical grades:
  mode=members in=grade1,grade2,grade3,grade4 inexact=t out=fresh_prefix
Writes ONLY .member_outcomes.tsv. Reads each .inputs.tsv, .overlap.raw.tsv and
.overlap.summary.tsv, including zero-call inputs; repeated input IDs across
partitions are rejected even when local query/ordinal numbers reset.
QUERY_EMITTED_MEMBER rows use the caller's member. REFERENCE_REPRESENTATIVE_MEMBER
rows use the reference's recorded member; a cross-member match remains on these
different axes. tp and tp_inexact are disjoint recorded matching classes;
matched is their sum. Nonzero endpoint offsets alone do not imply tp_inexact.
Unmatched/masked queries and missed/dont_care references remain separate, with
NA for inapplicable columns. Each axis partitions the union; do not sum axes or
interpret this as subtype accuracy. Reciprocal witnesses and totals are checked.
No legacy, timing or model-roster options are accepted in this summary-only mode.
Default mode=export remains the ordinary single-prefix exporter above.

Optional current-model winner distribution (canonical grades only):
  ref=main.fa,second.fa,third.fa expectedmodels=354,47,11 expectedinputs=10000
Adds .model_wins.tsv and .model_win_summary.tsv. Complete ordered FASTA rosters
retain models with zero wins. Each assigned scored QUERY match counts once for
its emitted model; inexact=t includes TP_INEXACT. Masks/unmatched calls do not
win, but every emitted model must still belong to the supplied roster.
Reports ALL and member min/Q1/median/Q3/max, zero wins, strictly fewer than15
and strictly fewer than46. Quantiles linearly interpolate sorted values at
(N-1)*P. These are winning matched calls on this panel, not construction-member
counts, unique model coverage, or proof that a zero-win model is dispensable.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.ExportCobalaminBenchmarkGrade "$@"

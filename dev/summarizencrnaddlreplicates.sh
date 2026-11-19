#!/bin/bash
# Report completed saved-window replicates with crossed seed/split uncertainty.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/summarizencrnaddlreplicates.sh in=workers.tsv out=fresh_prefix \
  seeds=12345,1,7,42 foldseeds=20261114,1,7 t=1 -Xmx512m -eoom

The complete worker roster has TAB-separated columns assignment_seed,ddl_seed,
table. Use assignment_seed=-1 for the fixed-not-seen check. Each table is the
out.replicates.tsv emitted by buildncrnaddllayoutcounts.sh study=replicates.
All49 methods and every judging fold must be present at every declared seed and
assignment. Only the first hashseed worker emits the unseeded index reference.
gc=t instead reports the48methods from study=replicategc, including freshly
measured w0 controls with the decided256/2048lists and60/100bonuses. Requires
the complete declared hashseed/assignment roster plus fixed check. Production
uses32hashseeds/tenassignments; smaller explicit rosters are finite witnesses.
Folds pool before nearest-rank p99, including rank509 all-acceptable-filtered
failures. No locus is dropped. Summary also reports wrongly removed acceptable
model/window pairs,windows losing all accepted models,missing-query-GC windows,
and loci losing every uniquely acceptable model. Window/pair counts are ALL only;
parent counts retain contribution strata. Fallback and zero-SD model/fold cases
are reported separately for the window and member statistics arms. Paired rows
compare each setting with its same-width/source w0,left-minus-right,with crossed SE.
Earlier followup/sizegrid/
definition rosters cannot be supplied with gc=t.
gccontrol=saved_size_worker.tsv in=new_gc_worker.tsv out=fresh_prefix instead
checks the exact1,024-bin256/2,048w0control against that prior matched worker.
Both source-arm controls must reproduce each original fold/population and all
six count/rank metrics. Writes out.parity.tsv. No full roster or new scoring.
Seeds and foldseeds are explicit,unique,nonnegative lists with at least2values.
Optional followup=workers.tsv adds the nearby-bonus/masked11 round. Its roster
columns are assignment_seed,ddl_seed,original_table,table. Every cell contains
the18new settings and two explicitly repeated100controls. Their exact per-fold
populations and six original metrics must match the paired original tables.
Exact508-bin overall/low rank histograms must reproduce count,sum,top10/50 and
contribution totals. Pool folds first,then use nearest-rank ceiling(0.99*N).
Reports add rank_p99 and low_rank_p99 for those20settings only;other old methods
have no saved distribution and are not rerun. Percentile uncertainty uses the
same crossed moments and paired cells as the other metrics.
Optional sizegrid=workers.tsv uses that same four-column schema and adds the
71size-job rows per cell;the matching256/4096winner is imported from followup.
Writes out.sizes.tsv with all72cells in both populations,rank/low/p99 means and
SDs,and differences from the best full-population cell with paired SEs. The
chosen full-population size has the fewest total entries within one paired SE
of the lowest full mean;equal total sizes prefer fewer short entries.
Optional definitions=workers.tsv shortsize=N longsize=N adds the three later
definition controls. Those sizes must match that deterministic grid choice.
All three explicit repeated controls,including0/0,retain exact old-fold parity.

Five disjoint folds are summed within each assignment;the fixed check sums its
two frozen judging halves. Means and sample standard deviations use those pooled
observations. Balanced crossed random-effects ANOVA separates hashseed,assignment
and seed-by-assignment interaction variance. Negative component estimates are
clipped at zero. Standard error of the mean is the square root of
hashVariance/seeds + splitVariance/assignments + interactionVariance/(both).
Paired standard errors apply this same calculation to cell-wise differences,
not to methods treated as independent. Fixed checks have seed variation only.

Writes out.summary.tsv and out.paired.tsv. Positive paired differences mean
left minus right;lower rank sums and higher top-k counts are better. Consensus
retention fractions average equally across folds and models within each size
third. NA retention for the index reference is not a measured zero. Full
consensus sequences and the saved acceptance oracle remain in all experiments;
this report makes no caller or independent-consensus-holdout claim.
gcall=t summarizes64all-data GC alternatives with the same crossed estimates.
Optional gcallfollowup=new_workers.tsv imports those64saved methods once from
in=original_workers.tsv and adds only12new window-GC alternatives. Both rosters
use assignment_seed,ddl_seed,table and the complete same crossed/fixed design.
Original and new worker method ranges are disjoint;missing or duplicate cells,
or an old control supplied as an extension,fail before reporting. Imported
controls and SD4/6retain every saved observation. New paired rows compare with
the original same-width window w0;combinations also compare with their range-only
gate. No original scoring/report command is rerun and no GC statistics are refit.
Use gcallcontrol=accepted_gc_table.tsv in=new_all_data_table.tsv for exact
no-GC parity at both widths/sources,all folds/populations and rank histograms.
Use gcpower=power.gc_power.tsv in=accepted.summary.tsv out=fresh_prefix for
a lookup-only candidate-cost table. This reads existing window-GC rankings and
the original population's retained-model counts;it runs no new scoring.
Surviving loci cost their saved best-window rank among retained models. Failed
loci exhaust the smallest retained candidate list across their saved windows.
This optimistic one-window estimate simulates neither caps nor patience and
is not a direct alignment count. Unranked thresholds have NA inferred costs.
Writes out.gc_cost.tsv with explicit no-filter baselines and count/percent savings.
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
CP="${DIR}current/"
if [ -f "${DIR}bbtools.jar" ];then CP="${DIR}bbtools.jar";fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.SummarizeNcrnaDdlReplicates "$@"

#!/usr/bin/env bash
# Saved endpoint/padding diagnostic; no caller, alignment or genome sequence read.
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/analyzecobalamininputpadding.sh in=baseline_grade,candidate_grade[,grade] \
 names=baseline,candidate[,candidate] traces=baseline.trace,candidate.trace[,trace] \
 profiles=baseline.cm_endpoints.tsv,candidate.cm_endpoints.tsv[,profile] \
 library=consensus.fa padding=candidate.padding.tsv[,padding.tsv] \
 bindings=panel.sources.tsv audit=resolved_alias_audit.tsv expectedinputs=10000 \
 out=fresh_prefix

Reads accepted normalized .witnesses.tsv and .genomes.tsv files, main-member
version-four CM traces and independent invocation profiles. Uses exact recorded
input/reference IDs; no rematching or proposal-dependent reference selection.
The consensus name, order and full sequence must match each padding resource.
There is one padding resource per candidate. The default baseline is fixed20/20;
basepadding=fixed_resource.tsv declares another sequence-bound uniform base.
minselectioncalls=1 and minselectiongroups=1 retain the old diagnostic defaults.
Use explicit larger minima for a selected resource; alias groups count once.
With one candidate, basepadding=... and exportselected=t also write
 .selected.CANDIDATE_NAME.tsv, an explicit base plus supported selected end pads.
Only the original P98_PLUS10 fit can supply exported overrides. Unchanged pads,
unsupported fits, ties and thin selection data retain the declared base.
For a final export, fitall=t pools every eligible saved pair across both group
halves without changing source-group identities. It requires exportselected=t.
Support and strict squared-error rules remain the same; endpoint estimates use
in_sample_n, and independent_judgment_common is zero. This does not add missing
genomes or teacher observations and does not run a caller or covariance search.

Reports:
 .coverage.tsv: each full witness population, pairwise and fixed all-arm cohorts.
 .model_ends.tsv: one row per candidate model/end, including unused models.
 .harms.tsv: pairwise previously-close ends made worse, retaining both models,
   actual decisions, original spans and configured/contig-clipped padding.
 .estimates.tsv: fixed-cohort judgment-half endpoint splices for fixed20,
   full candidate,3prime-only and selection-half per-model/per-end policies.
 .costs.tsv: uncertain nominal dense-work projections on actual searched requests,
   including unselected requests. Fixed/no-search requests are not invented work.

The alias audit assigns whole resolved groups by parity:0 chooses,1 judges.
A candidate model-end wins only with lower summed squared endpoint error on
the selection half; exact ties and thin support retain the base. All reports retain
the fixed recorded reference/path membership. Minor members must have unchanged
endpoints and are retained in the metrics. Same-model and model-switch rows are
separate. This policy split is not an independent model-training holdout.
The historical default baseline used a different state-context fit. Such changes
cannot be attributed to padding alone, even where configured/clipped flanks are
equal. For a new explicit-base comparison, bind compatible teacher, state-context
and all other resources before reading the saved endpoints.
The estimates therefore include unchanged/changed configured and clipped-window
strata; they do not claim an exact prediction of a new20-base-context caller.

Model/end better/worse counts overlap when different requests move in opposite
directions. Count each row once for unique model-ends, and use its better/worse
counts for request-end weighting. The clipped less/equal/more20 columns compare
the candidate flank with fixed20 on the SAME candidate original interval.
Configured padding can exceed20 while contig clipping erases that increase.

Optional expectedharm5=count[,count] reconciles accepted pairwise harm counts.
Optional cmcpu=seconds[,seconds] supplies observed full-candidate CM-search CPU,
EXCLUDING anchors. The projection scales it by CM length times triangular input
span on the SAME attempted requests. Sparse/retry work and changed context make
this an uncertain proxy, not a measured speedup. Real coupled caller cells must
measure both accuracy and cost. Export is opt-in; no policy adoption occurs.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.AnalyzeCobalaminInputPadding "$@"

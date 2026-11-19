#!/bin/bash
# Cross-fit compact partition weights from saved counts; no genome or alignment.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/optimizencrnapartitionweights.sh in=profiles512.tsv in2=profiles1024.tsv \
  panels=panels.tsv selection=selection.tsv ref=consensus.fa identities=frozen_prefix \
  ranks=accepted_partition.loci.tsv.gz out=fresh_prefix expectedmodels=508 \
  expectedgroups=3153 t=1 -Xmx2g -eoom

Loads existing member partitions. Reconstructs and checks the mixed consensus
quarter, and also tries consensus10-only without rebuilding member votes.
The raw-resource mode has three score terms, with one total for consensus;
it does not retain which k-mer length won a bin. Supply
counts=previous.counts.tsv.gz to reuse that immutable cache instead of sketching.
All consumed cache/resource/reference files must be bound by the execution recipe.

For arbitrary two-sketch layouts, replace in=/in2= with:
  layouts=producer.layouts.tsv counts=producer.counts.tsv.gz \
  control=legacy_member384_consensus128_both
Keep panels=, selection=, ref=, identities= and ranks= bound to the same original
roster. Validates the complete metadata/count interchange and exact
legacy raw control ranks before fitting historical raw layouts. Native DDL-only
layouts instead require matching native objects/count metadata, without a raw
control or reconstruction. Reads counts again one layout at a time;
never reconstructs keys in this mode. Output rows use layout_id/member_bins and
carry the validated metadata into out.layouts.tsv. Byte fields are payload,
not heap usage; no new counts file is written in this mode.
Optional fitlayouts=id1,id2 selects complete layout searches by exact ID; omitted
means all layouts. Full count/control validation still runs, including unselected
rows. Output remains in manifest order. Use a fresh output prefix for each subset.
Version-two mixed sketches have one match count each. Separate short/long
weights exist only where the key/bin partition encodes that distinction.
Version-one saved counts remain historical evidence; they cannot reconstruct
plain mixed-key matching from origin-filtered counts.
Layout mode also writes timing.tsv: actual score evaluations, fitting-thread CPU
and elapsed seconds for fitting/report enqueue only, excluding input parsing,
control validation and asynchronous writer CPU. Full-job time is measured separately.
Each selected layout emits start/completion diagnostics.

Search: randomdraws=5000 seeded log-uniform draws in [.1,10] by default;
set randomdraws=200 only to reproduce the earlier recorded search.
The declared reference term is fixed to one internally (member-long for
partitioned members; member-total for mixed members). A single pooled count
has no free weight and evaluates the unit point once. Otherwise use
two fit-only local refinements around the best five, with a 3^d stencil at
log steps log(2)/2 and log(2)/10. Every reported weight vector is normalized.
Objective: NOT_SEEN parent sum min(first acceptable-model rank,384); unavailable
ranks cost384. Both opposite judge halves must improve over unit weights before
an all-data REFIT is reported. A refit is not independently judged.
Weights w0..w3 follow metadata score_terms in active-prefix order, with unused
slots zero; raw-resource mode uses member-short/member-long/consensus-total.
Writes trials.tsv.gz, fits.tsv and decision.tsv under out=prefix. Writes a new
counts.tsv.gz only when counts= is omitted; cache replay keeps the original file.
Flatness ranges describe evaluated fit points, not confidence intervals.
Layout mode writes comparison.tsv, pareto.tsv and choices.tsv. The front uses
judged fitted cost only when BOTH halves improve, otherwise unit cost. Reports
both low-identity available-rank sums and capped sums separately. Choices give
the smallest stored member payload within 1% and 3% of the best permitted cost.
Merge disjoint completed batch comparisons with dev/summarizencrnalayoutfits.sh;
never use an all-data refit cost as an independently judged result.
No caller accuracy, adoption, new identities or alignment is implied.
Comparison rows also retain permitted overall/low-identity within10/50 counts,
available judging-parent count and uncapped available rank sum. Compare these
only against a baseline with the identical available population; unavailable
parents contribute384 each to the separate capped objective.
HELP
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
cd "$(dirname "$SCRIPT")/.."
DIR="$(pwd)/"
popd > /dev/null
CP="${DIR}current/"
if [ -f "${DIR}bbtools.jar" ]; then CP="${DIR}bbtools.jar"; fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.OptimizeNcrnaPartitionWeights "$@"

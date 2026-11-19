#!/usr/bin/env bash
# Compare complete saved paths without calling genomes or assigning accuracy.
set -eo pipefail
if [[ $# == 0 || $1 == -h || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/ncrnapathcalldelta.sh in1=OFF.gff in2=ON.gff \
 beforeconfig=OFF.caller.config afterconfig=ON.caller.config \
 sources=panel.sources.tsv genomes=OFF.genomes.tsv reference=OFF.references.tsv \
 expectedinputs=1000 out=fresh_directory -Xmx2g -eoom

Both GFFs require ordered input-ordinal boundaries. Exact flat caller configs
must bind identical input order and their respective saved output paths.
Use sources=the five-column panel source roster. For the frozen historical
panel only, legacyfilenames=t replaces sources=: each basename minus .fna.gz
must match one unique normalized OFF input ID; the complete bijection is checked.

Exact added/removed rows include endpoint, model and score changes. They are NOT
physical detection gains/losses. A unique mutual same-strand reciprocal-half
overlap among unmatched same-kind features gives a coordinate-change pair or,
when coordinates agree, a metadata/model-change pair. Ambiguous links remain
unpaired. Unlinked changes have no same-kind positive overlap in the other path.
Clan and emitted-member totals are separate; duplicate multiplicity is preserved.
Split/merge candidates count calls with multiple reciprocal overlap partners,
including unchanged anchors and duplicates; these are not inferred gene histories.

Reference coverage is same-input/contig/strand POSITIVE OVERLAP by a clan call.
It is separate from canonical one-to-one grader results. Scored and masked
references are retained separately. Changed CDS counts are observations only,
never CDS accuracy. Neighbour tables retain all overlaps and nearest flanks
on both strands, with ties, from both paths; these are not causal claims.
Whole-path overlap tables include unchanged calls: each unordered physical pair
with at least one clan call is counted once, including opposite-strand overlaps
and duplicate features. Affected clan calls, summed bases, maximum and mean
retain their separate denominators. These are ordinary reported-coordinate
overlaps; the deferred reporting-versus-selection excess is a different metric.

Outputs: feature counts, exact changed rows, overlap links, neighbour rows,
reference coverage and changed-feature/reference links, complete within-path
clan/feature overlaps and their totals, plus input bindings.
Incomplete or malformed streams fail without a completion record. Changed
call counts, unmatched locations or reference losses do not stop a comparison.
HELP
 exit 0
fi
DIR=$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=2g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$DIR/current/" prok.NcrnaPathCallDelta "$@"

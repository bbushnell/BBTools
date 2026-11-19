#!/bin/bash
# Merge disjoint completed layout-fit comparisons without re-running any search.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/summarizencrnalayoutfits.sh in=batch1.comparison.tsv,batch2.comparison.tsv \
  layouts=canonical.layouts.tsv out=fresh_prefix expectedlayouts=242 t=1 -Xmx1g -eoom

Requires disjoint completed comparisons with the same family, model/judge
population, random-draw budget and seed. Requires the exact canonical layout ID
set, geometry, model count, stored payload and query form; matching counts alone
are insufficient. Rejects duplicate, substituted or incomplete rosters.
Writes comparison.tsv, pareto.tsv and choices.tsv. Stored member bytes include
everything declared by the source metadata; runtime consensus and lookup indexes
are excluded from STORED resource size.

Pareto cost is opposite-half judged fitted cost only when both judge gains are
positive; otherwise unit weights remain the permitted choice. Attempted fitted
cost and low-identity available-rank/capped sums remain separately visible.
Choices select the smallest payload within 1% and 3% of the best permitted cost.
These are development comparisons, not untouched-test or caller results.
Matched-seed reporting: seedroster=roster.tsv layouts=canonical.layouts.tsv
reference=recipe_id out=prefix expectedlayouts=49 (omit in=). The roster has
ddl_seed and absolute comparison_table columns;each table needs its producer's
ddl_seed.txt sidecar. Requires every recipe on every matched seed,preserves
populations and the weight-search protocol,and reports mean/min/max plus paired
improvement against the reference. The comparison row's old seed column is the
weight-search seed,not the DDL hash seed. Writes seed_summary.tsv and four
objective-specific seed_choices.tsv rows;mean ordering alone is not adoption.
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
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.SummarizeNcrnaLayoutFits "$@"

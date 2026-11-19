#!/bin/bash
# Report final-consensus identity histogram quantiles and strict-cut bounds.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/summarize_ncrna_build_stats.sh in=stats1.tsv,stats2.tsv labels=first,second out=fresh.tsv expectedmembers=N cutoffs=0.64,0.70,0.75 t=1 -Xmx256m -eoom
Reads TrnaConsensusBuilder outbuildstats histograms, not sequences. Each histogram
must conserve expectedmembers. Reports nearest-rank quartiles of identities
rounded to 0.001 and bounds on the fraction strictly below each cutoff. A cutoff
must be a histogram-bin center; its own bin remains uncertain. No alignment or
caller execution. These are final-consensus identities only when the producing
builder retained final-consensus placement traces; verify that source provenance.
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
parseJavaArgs --xmx=256m --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.SummarizeNcrnaBuildStats "$@"

#!/bin/bash
# Global class2 conservation and raw-attempt curves from completed shards.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/summarize_ncrna_rejected_windows.sh in=finished_shards.tsv misses=full_misses.genes.tsv.gz population=PREFIX out=FRESH expectedlibraries=3 expectedcaptures=30 expectedrows=111455 expectedoriginals=111828 t=1 -Xmx2g -eoom
Manifest: library<TAB>capture<TAB>result_prefix. Every class2 original must have
one owning shard and all30stage/condition rows. Original coordinates and subset
flags are rejoined from accepted metadata. Curves and histograms are recomputed
from per-gene winners; per-window primary/edge work is summed once from summaries.
Own-model availability, eligibility and identity distributions remain explicit.
Identity bins are0.02wide, spans exact integer bases, and unavailable attempts
stay in each denominator. These are raw-attempt eligibility curves, NOT final
caller recall or an adopted cutoff. No DNA or alignment is read or executed.
Writes .curves.tsv, .histograms.tsv, .own_models.tsv and .work.tsv.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.SummarizeNcrnaRejectedWindows "$@"

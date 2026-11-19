#!/bin/bash
# Exclusive full-roster miss classes from saved completed baseline reports.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/summarize_ncrna_baseline_misses.sh in=grade_manifest.tsv population=prepared_prefix out=fresh_prefix expectedgenes=N expectedtruncated=N expectedincomplete=N t=1 -Xmx1g -eoom
Manifest header: library<TAB>grade_prefix<TAB>capture_manifest.
Paths must be absolute. Capture manifests use the baseline capture manifest schema;
an empty capture_id is allowed only for one standalone legacy capture.
Writes .classes.tsv and .genes.tsv.gz without DNA, alignment or grader replay.
Final matched calls take precedence, then unique pre-path acceptance, then actual
captured overlap versus no captured overlap. Scheduled gate rejections alone do
not count as captured DNA. Four exclusive classes conserve every original gene.
Native-truncated and incomplete-coverage subsets are reported independently.
Original IDs/frame/contig/coordinates/strand and completed summary counters are
checked. Actual accepted event and final call spans distinguish reciprocal-overlap
geometry failures from absent acceptance; class2 must not be called an identity
rejection. Rejected identity/HBM values are not present in these input schemas;
this tool does not fabricate their distribution from accepted-only fields.
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
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.SummarizeNcrnaBaselineMisses "$@"

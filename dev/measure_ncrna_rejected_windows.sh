#!/bin/bash
# Raw all-model saved-window attempts, without a genome reader or setting change.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/measure_ncrna_rejected_windows.sh in=misses.genes.tsv.gz library=LABEL population=PREFIX grade=PREFIX capture=PREFIX captureid=SHARD consensus=FASTA out=FRESH expectedmodels=N expectedclass2=N expectedgenes=N expectedincomplete=N expectedtruncated=N expectedshortboth=N t=16 -Xmx2g -eoom
Reads one saved capture shard and its associations to the complete class2 roster.
Each selected saved window is measured against every model, once per stage.
Production primary Quantum uses model/query and window/reference. Physical-edge
clipping uses original contig bounds; interior window ends never permit clipping.
Raw primary and edge stages are separate; rescue/trim/path are not replayed.
Curves are measured attempt eligibility at identity0.60/0.62/0.64, reciprocal>=0.5,
and span109/no-span. They are NOT predicted final recall or a selected setting.
The edge curve additionally requires genuine model clipping and the existing
edge identity floor. It measures all models at eligible physical edges, including
diagnostic attempts that the baseline may not have reached after an acceptance.
Histograms use identity bins[0,.02),[.02,.04),..., with1in a separate inclusive
last bin; span bins are exact integer bases. Absent positive overlap is explicit.
The short-both-ends subset is native cmfrom-1=44 AND cmlen-cmto=20; incomplete
uses the original saved flag. Subsets overlap and are not additive populations.
The best raw/own-model/six eligible attempts are retained separately per gene
and stage. Overlapping window work is counted once, never summed as per-gene CPU.
Writes .best.tsv.gz, .curves.tsv, .histograms.tsv and .summary.tsv. Bulk stays native.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.MeasureNcrnaRejectedWindows "$@"

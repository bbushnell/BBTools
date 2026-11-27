#!/usr/bin/env bash
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/reportncrnaseedsupport.sh in=ordered_seeds.fa seedset=fixed_prefix.fa \
 genes=original1.fa,original2.fa counts=100,200 expectedinitial=1975 \
 expectedkeys=2719 out=fresh_prefix t=1 -Xmx1g -ea -eoom

Count all sense-oriented 17-mer occurrences and distinct original gene records
supporting each selected key. Repeated words count as occurrences, once per gene
for copy support, and once toward distinct-hit coverage. Duplicate record IDs
do not collapse distinct gene copies. The fixed initial keys must be an unchanged
ordered prefix; all added keys must occur in the declared complete originals.

In recorded selection order, the first two distinct keys present in each gene
receive missing-slot credits. An added key with zero credits did not fill any
remaining two-hit slot when selected; this can expose batch overfill. Credits
may exceed one for a key that helps several genes. No key is selected or changed.

Writes .keys.tsv and .summary.tsv. Original pool counts and before/after coverage
are distinct from scored genomic-reference detection; no genomic scan is run.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError \
 -cp "${DIR}current/" prok.ReportNcrnaSeedSupport "$@"

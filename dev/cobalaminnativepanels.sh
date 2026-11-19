#!/bin/bash
# Metadata-only nested panels and raw member-qualified protection truth.
set -eo pipefail
if [[ $# == 0 || $1 == --help ]];then
 cat <<'HELP'
Usage: dev/cobalaminnativepanels.sh in=prepared/inputs.tsv manifest=selected/MANIFEST.tsv counts=recount/inputs.tsv loci=native_loci.tsv.gz sourceids=prepared/source_ids.tsv out=fresh_directory expectedgenomes=40686 expected=91072,2270,1786 rawexpected=91242,2270,1786 panels=1000,10000 t=1 -Xmx2g -eoom
Retains the supplied complete genome permutation exactly. Writes nested genome
lists, source/NEW-GFF bindings and six-column raw member protection truth, including
native truncations. Never opens genome DNA or adjudicates reference-GFF subtypes.
panels.tsv separates raw and construction-positive genome and locus counts;
per-member scored-positive counts are unavailable because that reference is a union.
seed_coverage.tsv retains every admitted original ID without sequence deduplication.
Use it with dev/ncrna_seed_coverage.sh fold=ALL representatives=f (no metadata=).
TRAIN in this compatibility ledger means all resource inputs; there is no holdout.
Outputs are fresh, with complete.tsv written last; partial directories are failures.
Optional chunkstart=10000 chunksize=1000 completedprefix=panel10000.sources.tsv
first verifies every skipped source against the accepted prefix, then also writes
disjoint suffix genome/source lists and chunks.tsv. Chunk caller ordinals restart
at zero; original input IDs are retained. No genome DNA is copied or read.
HELP
 exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [[ -h $SCRIPT ]];do
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.CobalaminNativePanels "$@"

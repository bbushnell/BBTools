#!/bin/bash
# Partition accepted source-frame metadata; no DNA or original-call rewrite.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/prepare_ncrna_baseline_shards.sh in=accepted.frames.tsv out=fresh_directory shards=10 expectedgenes=N expectedframes=N t=1 -Xmx1g -eoom
Balances saved available nucleotide counts by largest accepted assembly group.
Every group stays together. Global frame identities and raw ledger rows survive;
zero-gene and unavailable frames retain ownership. No source availability census,
DNA read, gene extraction, caller, alignment or resource construction is repeated.
Writes shard_N.frames.tsv, shard_N.reference.tsv and shard_N.inputs.config for
each shard, ownership.tsv and a final shards.tsv with whole-population counts.
Combine each input config with the same reviewed flat caller settings. The full
original gene ledger remains the canonical grading population. After capture,
dev/summarize_ncrna_roster_baseline.sh captures=manifest.tsv grades all shards
together and fits global GC moments without averaging partial standard deviations.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.PrepareNcrnaBaselineShards "$@"

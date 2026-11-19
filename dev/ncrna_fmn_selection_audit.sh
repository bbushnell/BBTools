#!/bin/bash
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_fmn_selection_audit.sh in=original_corpus.genomes.tsv
 catalog=earlier_genomes.tsv lineages=accepted_lineage_inputs.tsv.gz
 taxonomy=original_bac120_taxonomy.tsv.gz hits=original_standard_evidence.hits.tsv.gz
 out=fresh_prefix expectedgenomes=39979 expectedevaluation=11769
 expectedoverlap=10666 expectedlineages=4395 expectedhits=932681 t=1 -Xmx512m -eoom

Audits the literal earlier-PGAP intersection and its Bacteroidota/Chryseobacterium
counts. Original native table/line/target/strand/bounds and completed cut_ga
commands validate every counted FMN inclusion; legacy printed-score GA_pass is
ignored. Species overlap compares evaluation clean-core positive/negative groups
with discovery groups. No genome sequence, new search or split is performed.
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
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaFmnSelectionAudit "$@"

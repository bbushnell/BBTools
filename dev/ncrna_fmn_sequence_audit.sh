#!/bin/bash
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_fmn_sequence_audit.sh in=accepted_lineage_inputs.tsv.gz
 cores=saved_members.extract.tsv.gz metadata=original_GTDB_metadata.tsv.gz
 census=original_discovery_genomes.tsv previous=accepted_cross_selected.tsv
 out=fresh_prefix expectedqueries=152 expecteddiscovery=47 t=1 -Xmx1g -eoom

Selects the complete Chryseobacterium discovery population plus one each
Elizabethkingia and Empedobacter; extracts the declared evaluation native FMN
loci and exact31mer matches both strands. Best length is one contiguous exact
substring of one probe, not summed kmer coverage or a biological FMN call.
Compares actual original-file bytes/bases/records with the saved scan census;
retains lowercase/undefined/GC counts and original GTDB assembly metadata.
No CM search, training split, caller or resource change. Existing outputs fail.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaFmnSequenceAudit "$@"

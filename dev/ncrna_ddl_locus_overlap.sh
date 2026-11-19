#!/bin/bash
# Audit exact sequence and corrected source-locus presence among actual voters.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_ddl_locus_overlap.sh in=saved_rank_prefix out=fresh_prefix members=members.fa ledger=train.members.tsv assignments=assignments.tsv census=membership.models.tsv summary=membership.summary.tsv ref=consensus.fa genomes=corpus.genomes.tsv groups=groups.inputs.tsv.gz assemblies=assembly.inputs.tsv proxy=corpus.loci.tsv.gz recovered=recovered.ledger.tsv overlays=first.overlay.tsv.gz,second.overlay.tsv.gz t=1 -Xmx1g -eoom
Reads frozen construction inputs and saved parent/window/rank metadata and FASTA.
No genome, caller, new alignment, library generation or split change.
Only actual kept members contribute exact DNA or source-locus evidence.
Complete gene sequences come only from fully containing saved windows, with both
orientations handled explicitly. Partial or undefined queries retain UNKNOWN.
Recovered final cores and admitted current overlays replace raw annotation bounds.
Source-locus presence requires a resolved assembly, identical contig and strand,
and greater than half overlap of both query and corrected member spans.
SEEN means exact sequence OR same source locus; NOT_SEEN requires both absent.
This describes observed training overlap, not a newly independent test set.
Writes .loci.tsv, .labels.tsv and .summary.tsv. Labels feed the saved-rank reporter.
Defaults: expectedinputs=1000 expectedgenomes=39979 expectedtargets=2156
expectedcommon=2007 expectedmodels=508 expectedmethods=148. Fresh complete outputs
only; overwrite, append and sampling are rejected.
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
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaDdlLocusOverlap "$@"

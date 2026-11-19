#!/bin/bash
# Summarize saved construction membership, model-end classes and length spread.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]];then
	cat <<'HELP'
Usage: dev/summarize_ncrna_library_members.sh in=members.fa.gz originals=original_calls.tsv.gz consensus=library1.fa.gz,library2.fa.gz assignments=assignments1.tsv.gz,assignments2.tsv.gz labels=first,second out=fresh_prefix expectedmembers=N expectedoriginals=N t=1 -Xmx1g -eoom
Reads saved construction FASTA, emitted consensus FASTA, assignments and the full
original-gene ledger. No source genome or alignment is processed. Reports model
purity for the exact missing44left/20right class, separate original multiplicity,
and member/consensus length quantiles within and unweighted across models.
Each library writes PREFIX.LABEL.summary.tsv and PREFIX.LABEL.models.tsv.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "${DIR}current/" prok.SummarizeNcrnaLibraryMembers "$@"

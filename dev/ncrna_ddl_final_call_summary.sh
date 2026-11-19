#!/bin/bash
# Summarize proposal provenance in final saved CallGenes output.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_ddl_final_call_summary.sh in=saved_final.gff out=fresh_prefix mode=public t=1 -Xmx512m -eoom
mode=public requires absent proposal fields and reports unclassified/NA.
mode=hybrid requires native/sketch stream plus measured native shared counts.
mode=sketch requires sketch_only stream and native support NA.
Writes .calls.tsv (every RF00174 final RNA row with its original physical line,
model, proposal fields and geometry) and .summary.tsv (stream/support counts).
Duplicate predictions are retained. Other model namespaces are not targets.
These are final GFF rows, not unique loci, verifier accepts or graded TPs/FPs.
Original physical lines can join BindNcrnaCalls.line_map; bare contig names
cannot establish the input. No caller, alignment, sequence reads or regrading.
The arm contract must be declared; missing experimental fields fail loudly.
No overwrite, append, sampling or multithreaded readers.
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
parseJavaArgs --xmx=512m --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaDdlFinalCallSummary "$@"

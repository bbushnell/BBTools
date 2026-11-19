#!/bin/bash
# Export caller inputs and separate exact-sequence labels from saved metadata.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/ncrna_ddl_caller_panel.sh mode=capture in=rank.inputs.tsv audits=locus.loci.tsv out=fresh_prefix expectedinputs=1000 t=1 -Xmx512m -eoom
For the expanded population use mode=combined in=report.inputs.tsv and
audits=original.loci.tsv,additional.loci.tsv expectedinputs=3153.
Reads accepted saved metadata only; no sequence reads, calling or alignments.
Writes .config (one literal in= per selected input), .inputs.tsv (new/original
slot binding), .strata.tsv (all target parents, exact_sequence and source-locus
states kept separate), and .summary.tsv. Zero-target inputs and UNKNOWN labels
are retained. Declared population and per-input target counts must agree.
The combined input is the existing NcrnaSavedRankSummary .inputs.tsv, including
one selected representative per resolved noncontributing assembly group.
No overwrite, append, sampling or multithreaded metadata processing.
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
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaDdlCallerPanel "$@"

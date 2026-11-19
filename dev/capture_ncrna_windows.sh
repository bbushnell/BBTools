#!/bin/bash
# Captures native RF00174 candidate windows before shortlisting.
set -eo pipefail
if [[ $# == 0 || "$1" == --help ]]; then
	cat <<'HELP'
Usage: dev/capture_ncrna_windows.sh caller=caller.config reference=reference.tsv out=fresh_prefix t=8 -Xmx4g -eoom
Optional replaycontext=t saves separate .context.tsv/.context.fa.gz sidecars
with real refinement flanks for the supported ordinary RF00174 replay route.
Default candidate-only capture retains its existing route support and outputs.
The caller config selects original inputs and must explicitly set taxonomy=f.
The reference manifest binds input_id, absolute fasta and absolute gff paths.
RF00174 must be enabled; expectedmodels=508 checks the loaded model count.
Reports: .windows.tsv, .windows.fa.gz, .inputs.tsv and .summary.tsv.
results=t additionally requires selected-format DDL profiles and saves final
.calls.gff plus .results.tsv, one completed row per scheduled target window,
and .final_calls.tsv preserving each final call's original input identity even
when different files reuse contig names.
Results retain accepted identity/model/rank/score, actual primary/edge, rescue
and other alignment counts, and window GC/ACGT counts. Caching is disabled.
Pre-resolution verifier accepts and final call rows are separate populations;
this tool does not join them to the complete original-gene recall roster.
Coordinates are zero-based inclusive, both caller-oriented and original-forward.
Physical-edge retries reuse primary captured bytes and are not separate events.
Runs the configured ordinary caller once. Does not run an exhaustive oracle.
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
parseJavaArgs --xmx=2g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.CaptureNcrnaWindows "$@"

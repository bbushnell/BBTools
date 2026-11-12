#!/bin/bash
usage(){
cat <<'HELP'
Usage: dev/ncrnachangeevidence.sh in=sweep.changes.tsv inputs=sweep.inputs.tsv policies=policies.tsv calls=sweep family=cobalamin_ out=evidence.tsv

Extracts original stage and final-GFF rows around every paired changed reference.
Uses the native input-report order and explicit joint-policy roster to locate
the per-input/policy sidecars. Reads no FASTA, runs no caller and changes no grade.
Both strands and every final feature type are retained in the inspected region.
The baseline is policy-roster index 0. SEED rows attach to every requested locus
on that contig regardless of hit positions; they do not prove a seed at the locus.
Opposite-strand rows are included deliberately as neighboring evidence.

padding=0       Added genomic inspection context on either side of the reference.
                This is not a matching or causal-classification threshold.
models=f        Include MODEL observations only when explicitly enabled.
out=            Fresh TSV. Original rows are percent-framed in the final column;
                source path and physical line preserve their exact native origin.
                Converted coordinates are one-based genomic closed intervals.
                Seed rows have -1 coordinate fields and retain their raw hit list.
                Each empty locus/policy context gets NO_SAVED_EVIDENCE with both
                searched paths, source_line=0 and no invented original row.
                The completion marker counts inspected and empty contexts/loci.

All required sidecar files are checked before output creation. Missing files fail
loudly; they are not empty evidence. A run is complete only after successful close
and the NCRNA_CHANGE_EVIDENCE_COMPLETE diagnostic. Other read/parse failures can
leave incomplete output and never produce that completion marker.

WINDOW acceptance is before snapshot resolution and final path selection. Missing
observations of reused discovery are not rejection evidence. This tool extracts
evidence; it does not assign causal explanations or select a production policy.
One processing thread. Existing output files are protected; append unsupported.
HELP
}
if [[ $# == 0 || "$1" == -h || "$1" == --help ]]; then
	usage
	exit 0
fi
pushd . > /dev/null
SCRIPT="${BASH_SOURCE[0]}"
while [ -h "$SCRIPT" ]; do
	cd "$(dirname "$SCRIPT")"
	nextScript="$(readlink "$(basename "$SCRIPT")")"
	SCRIPT="$nextScript"
done
DIR="$(cd "$(dirname "$SCRIPT")/.." && pwd)/"
popd > /dev/null
if [ -f "${DIR}bbtools.jar" ]; then
	CP="${DIR}bbtools.jar"
else
	CP="${DIR}current/"
fi
source "${DIR}javasetup.sh"
parseJavaArgs --xmx=1g --xms=128m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -XX:+ExitOnOutOfMemoryError -cp "$CP" prok.NcrnaChangeEvidence "$@"

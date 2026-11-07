#!/bin/bash
# Summarizes existing matched-locus coordinates; never calls or rematches genes.
set -eo pipefail
usage(){
	cat <<'HELP'
Written by Brian Bushnell and Raiden

Usage: ncrnaendpointsummary.sh in=grade.matches.tsv out=summary family=RF00059 overlap=reciprocal_half

Reads the nine-column .matches.tsv report from NcrnaGffFamilyGrader or
NcrnaCalibrationGrade. Selects one exact family label. Coordinates must share
the report's contig namespace. Validates paired intervals, signed offsets,
the declared overlap rule and one-to-one truth/call locus use.

in=                 Complete existing matched-locus TSV, optionally compressed.
out=                Prefix for .summary.tsv, .histogram.tsv, .radii.tsv and
                    .quantiles.tsv. Parent directory must exist.
family=             Exact family label in the report; no inferred relabeling.
overlap=            Required declaration: reciprocal_half or half_shorter.
radii=0,1,2,3,5,10,20,50,100
                    Distinct nonnegative symmetric coordinate radii in bases.
expected=           Optional accepted number of matched loci.
exact5=             Optional accepted count of exact 5-prime ends.
exact3=             Optional accepted count of exact 3-prime ends.
both=               Optional accepted count of both ends exact.
ow=f                Existing output files are protected by default.

Positive errors point downstream in transcript orientation. Quantiles use nearest
rank. Joint radius counts require BOTH ends of the SAME locus to be within the
radius. These are coordinate-reachability diagnostics, not window containment,
neural-network accuracy, detection recall or biological-boundary confirmation.
No caller runs, rematching, resource changes or staging occur.
HELP
}
if [ $# -eq 0 ] || [ "$1" = -h ] || [ "$1" = --help ]; then
	usage
	exit 0
fi
SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
while [ -h "$SCRIPT" ]; do
	DIR="$(dirname "$SCRIPT")"
	RESOLVED_SCRIPT="$(readlink "$SCRIPT")"
	SCRIPT="$RESOLVED_SCRIPT"
	if [ "${SCRIPT#/}" = "$SCRIPT" ]; then
		RESOLVED_SCRIPT="$DIR/$SCRIPT"
		SCRIPT="$RESOLVED_SCRIPT"
	fi
done
DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
if [ -f "$DIR/bbtools.jar" ]; then
	CP="$DIR/bbtools.jar"
else
	CP="$DIR/current/"
fi
source "$DIR/javasetup.sh"
parseJavaArgs --xmx=1g --xms=64m --mode=fixed "$@"
setEnvironment
exec java $EA $EOOM $SIMD $XMX $XMS $XSS -cp "$CP" prok.NcrnaEndpointSummary "$@"

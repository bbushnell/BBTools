#!/bin/bash

usage(){
echo "
Written by Eru
Last modified July 31, 2026

Description:  Computes the CheckM1-style marker-counting completeness and
contamination estimate for a genome bin from its marker vector.  This is the
ORACLE/baseline the MAG-QC completeness/contamination neural net is validated
against (the shipped estimator is the net; this is the counting method it is
checked against).

Completeness = detected marker families / expected marker families.
Contamination (headline) = excess copies / expected families, where excess
copies = sum over families of max(0, copies-1) (so 3 copies counts more than 2).
A simpler multi-copy-family contamination is also reported.

Usage:  magqc.sh vector=<vec.tsv> out=<report.tsv>

The vector file is the marker-vector TSV written by 'markervector.sh out=<file>'
(family_id / representative / count rows plus '#'-prefixed derived scalars).

Required parameters:
vector=<file>   Marker-vector TSV from markervector.sh.

Optional parameters (and their defaults):
out=stdout      Output report TSV.
ow=t            (overwrite) Overwrite existing output.

Output ('metric<tab>value' lines):
completeness_pct, contamination_pct (headline, excess-copy),
contamination_pct_multicopy (secondary), expected_markers, detected_markers,
multicopy_markers, excess_copies, effective_denominator, domain_assignment,
marker_set_id, lineage_taxid, rank, assignment_confidence, ood_status
(fixed 'unknown' until OD-9), sufficient_evidence.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
-eoom           Exit if an out-of-memory exception occurs.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"

}

if [ -z "$1" ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
	usage
	exit
fi

resolveSymlinks(){
	SCRIPT="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	if [ -f "$DIR/bbtools.jar" ]; then
		CP="$DIR/bbtools.jar"
	else
		CP="$DIR/current/"
	fi
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=2g" "--xms=2g" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prot.MagQCCLI $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"

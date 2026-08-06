#!/bin/bash

usage(){
echo "
Written by Eru
Last modified July 31, 2026

Description:  Turns a genome bin's proteins into a fixed-length gene
presence/count vector against a marker set.  For each selected single-copy
marker family, entry i = how many of the bin's proteins match that family
(the copy count; 0 = absent).  Each bin protein is aligned (BLOSUM62) against
the marker family representatives and assigned to its single best-matching
family above the identity/coverage thresholds.  This is the feature-extraction
step that feeds the MAG-QC completeness/contamination neural net; it does not
run the net itself.

Usage:  markervector.sh bin=<bin.faa> markers=<markers.faa> out=<vec.tsv>

The markers file is the marker-representatives FASTA written by
'markerfactory.sh repsout=<file>' (self-describing headers carrying family_id,
domain, selected, genomes and copies).  The plain marker TSV lacks the
representative sequences and cannot be scored against on its own.

Required parameters:
bin=<file>      Genome bin protein FASTA (amino acids).
markers=<file>  Marker-representatives FASTA (markerfactory repsout=).

Optional parameters (and their defaults):
out=stdout      Output TSV.
domain=         Which domain's marker set to score against (required only if the
                marker file holds more than one domain).
minid=          Override min percent identity to assign a protein to a family
                (default: the marker set's build identity, else 90; 0.9 also ok).
mincov=         Override min aligned fraction of protein and rep
                (default: the marker set's build coverage, else 0.8).
ow=t            (overwrite) Overwrite existing output.

Output columns (tab-separated):
family_id representative count
One row per selected marker family, followed by '#'-prefixed derived scalars
(families present / exactly-once / multi-copy, proteins matched / unmatched).

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prot.MarkerVectorCLI $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"

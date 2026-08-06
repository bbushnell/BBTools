#!/bin/bash

usage(){
echo "
Written by Eru
Last modified July 31, 2026

Description:  Greedy identity-threshold protein clustering (CD-HIT / linclust
shape).  Reads a protein FASTA and writes a representative-to-member TSV, one
row per member.  Sequences are processed longest-first; each is seeded by
shared k-mers against existing representatives, aligned with BLOSUM62
affine-gap local alignment (gap-open 11, gap-extend 1), and joins the best
representative meeting the identity and coverage thresholds, or starts a new
cluster.

Usage:  clusterproteins.sh in=<proteins.faa> out=<clusters.tsv>

Required parameters:
in=<file>       Protein FASTA (amino acids).

Optional parameters (and their defaults):
out=stdout      Output TSV.  A .meta sidecar is written beside it.
minid=90        Minimum percent identity to join a cluster (0.9 also accepted).
mincov=0.8      Minimum aligned fraction of both member and representative.
k=5             Seed k-mer length.
reduced=f       Use amino8 reduced-alphabet seeds (more sensitive).
ow=t            (overwrite) Overwrite existing output.

Output columns (tab-separated):
cluster_id representative member identity coverage is_representative

Notes:
One row per member; the representative is listed as a member of its own
cluster.  Cluster ids are stable within a run only (cross-run stable ids and
incremental update are not implemented).

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
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP prot.ClusterProteins $@"
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"

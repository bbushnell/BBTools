#!/bin/bash

usage(){
echo "
Written by Noire and Brian Bushnell
Last modified July 26, 2026

Description:  Builds neural-network training vectors for QuickClade hit RANKING from a
QuickClade machine-format hit table.  One 33-dimension vector per hit plus a regression
label; the trained network reorders a query's hits, replacing the compositeScore heuristic.

The input MUST come from a run with printcomposite=t, which appends the trailing
CompositeScore column; the tool exits with an error rather than scoring a table without it.
Vectors are emitted for the top maxemit hits of each query.

Label is the LCA depth of the hit against the query's true taxID, on an 11-rung ladder:
same_node=1.0, species=0.9, genus=0.8 ... domain=0.1, no agreement=0.  The query's true
taxID is read from the shred header (tid_NNN); a query lacking one is skipped rather than
given a fabricated label.

Usage:  rankingvectorizer.sh in=hits.tsv out=vectors.tsv tree=tree.taxtree.gz buckets=4096

Generating suitable input:
quickclade.sh <shreds> sketchfile=4k server=f slow callssu=t records=10 heapsize=10 \\
  printcomposite=t percontig format=machine usetree out=hits.tsv

Parameters:
in=<file>       QuickClade machine-format hit table, produced with printcomposite=t.
out=<file>      Output vector file, with a tab-delimited #dims header.
tree=<file>     Taxonomy tree (tree.taxtree.gz); required, for LCA and lineage.
buckets=0       DDL sketch bucket count of the reference DB used for the run: 0 for fast
                mode (no sketch), 4096 for sketchfile=4k, 32768 for sketchfile=32k.
                This feeds an input dimension, so a wrong value trains on a wrong feature.
maxemit=10      Emit vectors for this many top hits per query.  Governs row alignment for
                every downstream tool, so keep it consistent across a training set.
eps=0.02        Stabilizer in the normalized-composite dimension, (x-t)/((m-t)+eps).
                Measured, not arbitrary: below a within-query spread of 0.02 the composite
                ordering is at or below a coin flip in every clade tested, so without eps
                the encoding stretches a pure noise band across the full output range and
                invents an ordering.  Do not lower it without repeating that measurement.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
-da             Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
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

	parseJavaArgs "--xmx=8g" "--xms=8g" "--percent=42" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP clade.RankingVectorizer $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"

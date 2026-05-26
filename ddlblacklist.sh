#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Noire
Last modified May 24, 2026

Description:  Builds a DDL kmer blacklist from pre-built DDL sketch files.
Reads sketches with kmer arrays (built with ddlwriter.sh kmers=t), counts
how many distinct records each kmer appears in, and outputs kmers exceeding
the threshold as FASTA.  Low-entropy kmers can also be blacklisted.

Usage:  ddlblacklist.sh in=sketches.tsv out=blacklist.fa [mintaxcount=20]

Parameters:
in=<file>       Input DDL sketch file (must have been built with kmers=t).
out=<file>      Output FASTA file of blacklisted kmer sequences.
mintaxcount=20  Minimum number of records (taxa) a kmer must appear in
                to be blacklisted.
minentropy=0    Minimum Shannon entropy threshold (0-1).  Kmers below this
                entropy are blacklisted regardless of taxa count.
                BBSketch default is 0.66 with entropyk=3.
entropyk=3      Sub-kmer length for entropy calculation.
k=19            K-mer length (must match the input sketch file).

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will
                specify 200 megs. The max is typically 85% of physical memory.
-eoom           This flag will cause the process to exit if an out-of-memory
                exception occurs.  Requires Java 8u92+.
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

	parseJavaArgs "--xmx=4g" "--xms=4g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ddl.DDLBlacklistMaker $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"

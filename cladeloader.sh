#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified October 13, 2025

Description:  Loads fasta files and writes clade files.

Usage: cladeloader.sh in=contigs.fa out=clades.clade

Parameters:
in=<file,file>  Fasta files with tid in headers.
out=<file>      Output file.
maxk=5          Limit max kmer length (range 3-5).
a48             Output counts in ASCII-48 instead of decimal.
16s=<file,file> Optional tax-labeled file of 16S sequences.
18s=<file,file> Optional tax-labeled file of 16S sequences.
replaceribo     Set true if existing ssu should be replaced by new ones.
usetree=f       Load a taxonomic tree to generate lineage strings.
aligner=quantum Options include ssa2, glocal, drifting, banded, crosscut.

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
	SCRIPT="$0"
	while [ -h "$SCRIPT" ]; do
		DIR="$(dirname "$SCRIPT")"
		SCRIPT="$(readlink "$SCRIPT")"
		[ "${SCRIPT#/}" = "$SCRIPT" ] && SCRIPT="$DIR/$SCRIPT"
	done
	DIR="$(cd "$(dirname "$SCRIPT")" && pwd)"
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=4g" "--xms=4g" "--percent=42" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP clade.CladeLoader $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"

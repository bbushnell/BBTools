#!/bin/bash

usage(){
echo "
Written by Brian Bushnell and Noire
Last modified May 19, 2026

Description:  Identifies and classifies ribosomal SSU (16S/18S) sequences
using DynamicDemiLog (DDL) sketching against a pre-built SSU reference
database.  Two modes:
  Default: input sequences are assumed to be SSUs.  Each is classified as
    16S or 18S by alignment against consensus sequences, sketched, and
    compared to the reference database.
  Call mode: input is genomic sequence.  Gene-calling finds all SSUs
    (potentially multiple per contig), then classifies and compares each.

Usage:  findssu.sh ssu1.fa [ssu2.fa ...] ref=<ssu_ddls.tsv> [records=5]
    or: findssu.sh genome.fa call ref=<ssu_ddls.tsv>
    or: findssu.sh ssu.fa ref16s=<16S.tsv> ref18s=<18S.tsv>

SSU Mode (default):
Each input sequence is aligned to 16S and 18S consensus sequences to
determine type, then sketched with DDL and compared to the reference.

Call Mode:
Gene-calling identifies all SSU sequences in the input genome(s).
Each SSU found is individually sketched and classified.

Parameters:
ref=<file>      Pre-built SSU DDL reference file (TSV format).
ref16s=<file>   Separate 16S reference file.
ref18s=<file>   Separate 18S reference file.
call            Enable gene-calling mode for genomic input.
records=5       Max hits to display per query.
minhits=1       Minimum matching DDL buckets to report a hit.
index=f         Use inverted index for query acceleration.
k=13            K-mer length for hashing.
buckets=128     Number of DDL buckets.
exponent=5      Exponent bits.
t=1             Number of threads.

Output columns:
ANI, WKID, Complt, Matches, Type, qLen, rLen, TID, Query, Name,
File, Contig, Start, Strand

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

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP ddl.SSUCompare $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"

#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified December 9, 2025

Description:  Fast lightweight scanner that parses sequence files.
Reports base and record counts.  Performs basic integrity checks;
reports corruption and exits with code 1 when detected.

Usage:  fastqscan.sh <file> <threads>
e.g.
fastqscan.sh contigs.fasta
fastqscan.sh reads.fq.gz
fastqscan.sh reads.fq 2

Input may be fastq, fasta, sam, scarf, gfa, or fastg, 
compressed or uncompressed.  To input stdin use e.g. stdin.fq
as the argument (with proper extension).
FastqScan does not perform validation of all fields with the same rigor as
Reformat or other BBTools, but it does detect typical problems like incomplete
records, missing a + or @ symbol, missing fields, and base/quality length
mismatches, and reports Windows-style CRLF newlines.  It can perform 
multithreaded reading for fastq only, which disables integrity checks.
Bgzipped input processing is multithreaded and far faster than regular gzip.
SIMD support is autodetected and can be disabled with the flag simd=f.

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
	CP="$DIR/current/"
}

EA="-ea"
SIMD="--add-modules jdk.incubator.vector"
XMX="-Xmx1g"
XMS="-Xms128m"

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=1g" "--xms=128m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP stream.FastqScan $@"
	#echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
#setEnv "$@"
launch "$@"

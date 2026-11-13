#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 14, 2024

Description: Counts and summarizes the number of reads with each barcode,
using class BarcodeStats.  Can also do barcode assignment.

Usage:   countbarcodes2.sh in=<file> counts=<file>

Input may be stdin or a fasta or fastq file, raw or gzipped.

Input Parameters:
in=<file>           Input reads, whose names end in a colon then barcode.
countsin=<file>     Input of counts; optional.
quantset=<file>     Only quantify barcodes in this file.
interleaved=auto    (int) If true, fastq input will be considered interleaved.
expected=           Comma-delimited list of expected bar codes.

Output parameters:
maxrows=-1          Optionally limit the number of rows printed.
printheader=t       Print a header.
out=<file>          (counts) Write bar codes and counts here.  'out=stdout' 
                    will pipe to standard out.
barcodesout=<file>  Barcode assignment counts.
mapout=<file>       Map of observed to expected barcode assignments.
outcontam=<file>    Requires labeled data, and causes contam quantification.

Processing Parameters:
countundefined=t    Count barcodes that contain non-ACGT symbols.
pcrmatrix=f         Use a PCRMatrix for barcode assignment.
mode=hdist          PCRMatrix type.
mincountpercentile=0
                    Fraction [0,1]. If positive, keep only observed barcodes
                    at or above this count-mass-weighted percentile cutoff for
                    matrix/server assignment; ties at the cutoff are retained.
                    For example, counts 1,1,8,90 with mincountpercentile=0.1
                    keep 8 and 90. Default 0 disables this filter. Raw counts
                    and raw statistics are unchanged whether this is enabled or
                    disabled.


Java Parameters:
-Xmx                This will set Java's memory usage, overriding autodetection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
-eoom               This flag will cause the process to exit if an
                    out-of-memory exception occurs.  Requires Java 8u92+.
-da                 Disable assertions.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
For documentation and the latest version, visit: https://bbmap.org
"
}

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
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

	parseJavaArgs "--xmx=2g" "--xms=2g" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP barcode.CountBarcodes2 $@"
	echo "$CMD" >&2
	java $EA $EOOM $SIMD $XMX $XMS -cp "$CP" barcode.CountBarcodes2 "$@"
}

resolveSymlinks
setEnv "$@"
launch "$@"

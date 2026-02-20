#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 16, 2024

Description:  Accepts multiple input files from a demultiplexed lane.
Parses the barcode from the filename and adds (tab)BARCODE to read headers.
Outputs all reads into a single file.  Optionally, trims bases and drops R2.
Intended for evaluating demultiplexing methods.  For example:
tagandmerge.sh path/*0.*.fastq.gz dropr2 trim out=tagged.fq.gz barcodes=bc.txt

Usage:  tagandmerge.sh *.fastq.gz out=<output file>
or
tagandmerge.sh in=<file,file,file> out=<output file>

Input may be fasta or fastq, compressed or uncompressed.

Standard parameters:
in=<file,file>  A comma-delimited list of files.  If wildcards are used,
                omit in= and the commas.
out=<file>      Print all reads to this destination.
barcodes=<file> Print barcodes from file names to this destination.
trim=-1         If positive, trim all reads to this length.
dropr2=f        Discard read 2 if the input is interleaved.
shrinkheader=f  (shrink) Illumina only; remove unnecessary header fields.
remap=-+        Remap symbols in the barcode.  By default, '+' replaces '-'.
                To eliminate this set 'remap=null'.

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

	parseJavaArgs "--xmx=300m" "--xms=300m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP barcode.TagAndMerge $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"

#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 17, 2015

Description:  Concatenates barcodes and quality onto read names.

Usage:        mergebarcodes.sh in=<file> out=<file> barcode=<file>

Input may be stdin or a fasta or fastq file, raw or gzipped.
If you pipe via stdin/stdout, please include the file type; e.g. for gzipped fasta input, set in=stdin.fa.gz

Optional parameters (and their defaults)

Input parameters:
in=<file>       Input reads. 'in=stdin.fq' will pipe from standard in.
bar=<file>      File containing barcodes.
int=auto        (interleaved) If true, forces fastq input to be paired and interleaved.
qin=auto        ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.

Output parameters:
out=<file>      Write muxed sequences here.  'out=stdout.fa' will pipe to standard out.
overwrite=t     (ow) Set to false to force the program to abort rather than overwrite an existing file.
ziplevel=2      (zl) Set to 1 (lowest) through 9 (max) to change compression level; lower compression is faster.
qout=auto       ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).

Other parameters:
pigz=t          Use pigz to compress.  If argument is a number, that will set the number of pigz threads.
unpigz=t        Use pigz to decompress.

Java Parameters:
-Xmx            This will set Java's memory usage, overriding autodetection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.
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
	CP="$DIR/current/"
}

setEnv(){
	. "$DIR/javasetup.sh"
	. "$DIR/memdetect.sh"

	parseJavaArgs "--xmx=3200m" "--xms=3200m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.MergeBarcodes $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"

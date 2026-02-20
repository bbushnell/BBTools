#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 28, 2026

Description:  Converts between sam, bam, fasta, fastq.
              Supports subsampling, paired files, and multithreading.

Usage:  stream.sh in=<file> out=<file> <other arguments>
or
stream.sh <input_file> <output_file> <other arguments>
e.g.
stream.sh mapped.bam mapped.sam.gz
stream.sh in=reads.fq out=subset.fq samplerate=0.1

File parameters:
in=<file>       Primary input file, type detected from extension.
in2=<file>      Secondary input file for paired reads.
out=<file>      Primary output file, optional, type based on extension.
out2=<file>     Secondary output file for paired reads.
                Note: Use # symbol for auto-numbering, e.g. reads_#.fq

Processing parameters:
samplerate=1.0  Fraction of reads to keep (0.0 to 1.0).
sampleseed=17   Random seed for subsampling (-1 for random).
reads=-1        Quit after processing this many reads (-1 = all).
ordered=t       Maintain input order in output.
skipreads=0     Skip this many initial reads (or pairs).

Threading parameters:
threadsin=-1    Reader threads (-1 = auto).
threadsout=-1   Writer threads (-1 = auto).

Other parameters:
simd            Add this flag for turbo speed. Requires Java 17+ and AVX2,
                or other 256-bit vector instruction sets.

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

	parseJavaArgs "--xmx=2g" "--xms=2g" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP stream.StreamerWrapper $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
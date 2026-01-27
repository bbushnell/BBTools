#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified January 21, 2026

Description:  Makes cov files for QuickBin.
Notes: This program can use a lot of memory when there are many contigs.
Sam files use less memory than bam.  

Usage:
covmaker.sh *.sam out=cov.txt
or
covmaker.sh in=sample1.bam,sample2.bam out=cov.txt
or
covmaker.sh cov.txt out=cov7.txt condense=7

File parameters:
in=<file>       Input.  Properly named files (*.sam, etc) do not need 'in='.
                Sam, bam, and cov files are supported as input;
		bam uses the most memory, and cov the least.
condense=       When there are more than this many samples (sam/bam files),
                combine some into the same logical sample to save memory.
reorder=t       Reorder samples by decreasing entropy to improve indexing.
readthreads=4   Load up to this many sam/bam files concurrently.
                Lower uses less memory (when there are more samples).
mincontig=100   Ignore contigs shorter than this.  Saves memory.
minseed=2.5k    Don't calculate entropy from contigs shorter than this.

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

	parseJavaArgs "--xmx=4000m" "--xms=4000m" "--percent=84" "--mode=auto" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP bin.CovMaker $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"

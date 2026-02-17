#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified February 16, 2026

Description:  Makes cov files for QuickBin.
Notes: This program can use a lot of memory when there are many contigs.
Sam files use less memory than bam.  

Usage:
covmaker.sh *.sam out=cov.txt
or
covmaker.sh in=sample1.bam,sample2.bam out=cov.txt
or
covmaker.sh cov.txt out=cov7.txt condense=7 reorder

File parameters:
in=<file>       Input.  Properly named files (*.sam, etc) do not need 'in='.
                Sam, bam, and cov files are supported as input;
		bam uses the most memory, and cov the least.
out=<file>      Output coverage file.

Other parameters
condense=<int>  When there are more than this many samples (sam/bam files),
                combine some into the same logical sample to save memory.
reorder=t       Reorder samples by decreasing entropy to improve indexing.
mincontig=100   Ignore contigs shorter than this.  Saves memory.
readthreads=4   Load up to this many sam/bam files concurrently.
                Lower uses less memory (when there are more samples).
minseed=2.5k    Don't calculate depth entropy from contigs shorter than this.
magnitude=t     Use magnitude (total sample volume) in merge decisions -
                prioritize low-volume samples.
cosine=t        Use cosine similarity in merge decisions - 
                prioritize similar samples.
entropy=f       Use depth entropy in merge decisions - 
                prioritize low-entropy samples.
negcos=f        Invert the cosine function to prioritize dissimilar samples.
magpower=1.0    Raise sample volume to this power to alter its strength.
entpower=1.0    Raise entropy to this power.
compare=100k    Only compare this many largest contigs when calculating 
                pairwise depth similarity.
lognorm=f       Use logs of normalized depth, instead of raw depth,
                for pairwise similarity.

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

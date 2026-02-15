#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified August 22, 2024

Description:  Filters a sam file to select only reads with substitution errors
for bases with quality scores in a certain interval.  Used for manually
examining specific reads that may have incorrectly calibrated quality scores.

Usage:  filtersubs.sh in=<file> out=<file> minq=<number> maxq=<number>

Parameters:
in=<file>       Input sam or bam file.
out=<file>      Output file.
minq=0          Keep only reads with substitutions of at least this quality.
maxq=99         Keep only reads with substitutions of at most this quality.
countindels=t   Also keep reads with indels in the quality range.
minsubs=1       Require at least this many substitutions.
minclips=0      Discard reads with more clip operations than this.
maxclips=-1     If nonnegative, discard reads with more clip operations.
keepperfect=f   Also keep error-free reads.

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

	parseJavaArgs "--xmx=120m" "--xms=120m" "--mode=fixed" "$@"
	setEnvironment
}

launch() {
	CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP jgi.FilterReadsWithSubs $@"
	echo "$CMD" >&2
	eval $CMD
}

resolveSymlinks
setEnv "$@"
launch "$@"
